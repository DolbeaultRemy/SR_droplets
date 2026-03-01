function objective = evaluateFeatureDetectionPipeline(x, img, gtMask, params, doPlot)
%% evaluateFeatureDetectionPipeline
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Objective wrapper for Bayesian optimization.
%   Runs the detection + shape extraction pipeline for a candidate
%   parameter set and returns a score relative to ground truth.
%
% Inputs:
%   x        - struct or table from bayesopt with tunable params
%   img      - 2D image (double or uint)
%   gtMask   - logical mask of ground truth region(s)
%   params   - base parameter struct (fixed params)
%   doPlot   - optional boolean flag to display results (default=false)
%
% Output:
%   objective - scalar (to be MINIMIZED by bayesopt)
%
% Notes:
%   Optional notes, references.

    if nargin < 5
        doPlot = false;
    end

    % --- Copy params and overwrite tunable fields ---
    p = params;

    % Clip / enforce constraints
    p.dogGaussianSmallSigma    = min(max(x.dogGaussianSmallSigma,0.01), x.dogGaussianLargeSigma-0.01);
    p.dogGaussianLargeSigma    = max(x.dogGaussianLargeSigma, p.dogGaussianSmallSigma+0.01);
    p.adaptiveSensitivity      = min(max(x.adaptiveSensitivity,0),1);

    % Use nearestOdd helper for neighborhood size
    p.adaptiveNeighborhoodSize = nearestOdd(x.adaptiveNeighborhoodSize, 3);

    p.minPeakFraction          = min(max(x.minPeakFraction,0),1);
    p.minimumPatchArea         = max(round(x.minimumPatchArea),1);
    p.intensityThreshFraction  = min(max(x.intensityThreshFraction,0),1);
    p.edgeSigma                = max(x.edgeSigma,0.01);

    lowT  = min(max(x.edgeThresholdLow,0),1);
    highT = min(max(x.edgeThresholdHigh,0),1);
    if highT <= lowT, highT = min(lowT+0.05,1); end
    p.edgeThresholdLow  = lowT;
    p.edgeThresholdHigh = highT;

    % --- Run detectPatches safely ---
    try
        [patchProps, ~, imgCropped, xStart, yStart] = Analyzer.detectPatches(img, p);
    catch ME
        warning('detectPatches failed: %s', ME.message);
        objective = 0;
        return;
    end

    % --- Run extractAndClassifyShapes safely ---
    try
        results = Analyzer.extractAndClassifyShapes(imgCropped, patchProps, xStart, yStart, p);
    catch ME
        warning('extractAndClassifyShapes failed: %s', ME.message);
        objective = 0;
        return;
    end

    % --- Build predicted mask ---
    predictedMask = false(size(img));
    for k = 1:numel(results)
        BW = results(k).BW;
        if isempty(BW), continue; end
        [h,w] = size(BW);
        yIdx = (yStart:yStart+h-1);
        xIdx = (xStart:xStart+w-1);
        predictedMask(yIdx, xIdx) = predictedMask(yIdx, xIdx) | BW;
    end

    % --- Compute Dice score ---
    if any(gtMask(:))
        score = dice(predictedMask, gtMask);
    else
        score = 0;
    end

    % --- Optional plotting ---
    if doPlot
        p.gtMask = gtMask;  % overlay
        Plotter.plotDetectedPatches(img, patchProps, results, p, xStart, yStart, 'Optimization Candidate');
        drawnow;
    end

    % bayesopt minimizes → return negative Dice
    objective = -score;
end

function val = nearestOdd(x, minVal)
%NEARESTODD  Round to the nearest odd integer, enforcing a minimum.
%
%   val = NEARESTODD(x) rounds x to the nearest odd integer.
%   val = NEARESTODD(x, minVal) enforces that the result is at least minVal.

    if nargin < 2
        minVal = 1;
    end

    % Round to nearest integer
    n = round(x);

    if mod(n,2) == 0
        lowerOdd = n - 1;
        upperOdd = n + 1;

        % Choose whichever odd is closer to original x
        if abs(lowerOdd - x) <= abs(upperOdd - x)
            n = lowerOdd;
        else
            n = upperOdd;
        end
    end

    % Enforce minimum
    val = max(n, minVal);
    if mod(val,2) == 0
        val = val + 1; % bump to next odd
    end
end
