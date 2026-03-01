function results = analyzeAutocorrelation(autocorrresults)
%% analyzeAutocorrelation
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Computes features from g2(theta) curves and aggregates results.
%
% Notes:
%   Assumes theta_values are given in radians.

    feature_list     = table();
    per_group        = cell(numel(autocorrresults.scan_parameter_values),1);
    per_group_mean   = cell(numel(autocorrresults.scan_parameter_values),1); % mean-curve features

    % --- Loop over groups ---
    for i = 1:numel(autocorrresults.scan_parameter_values)
        curves  = autocorrresults.g2_curves{i};
        thetas  = autocorrresults.theta_values;

        group_tbl = table();
        for j = 1:size(curves,1)
            g2curve      = curves(j,:);
            feats        = extractFeatures(thetas, g2curve);
            group_tbl    = [group_tbl; feats]; %#ok<AGROW>
            feature_list = [feature_list; feats]; %#ok<AGROW>
        end
        per_group{i} = group_tbl;

        % --- Compute mean-curve features for this group ---
        meanFeatures = extractMeanG2Features(thetas, curves);
        per_group_mean{i} = meanFeatures;
    end

    %% Compute group summaries (mean & SEM of features)
    summary_table = table();
    N_params = numel(autocorrresults.scan_parameter_values);

    for i = 1:N_params
        FT             = per_group{i};

        % Scalars
        vars           = {'A2','A4','A6','S2','Q4','H6','Contrast'};
        
        % Group averages
        meanVals       = varfun(@mean, FT, 'InputVariables', vars);
        semVals        = varfun(@(x) std(x,0,1,'omitnan')/sqrt(size(FT,1)), ...
                          FT, 'InputVariables', vars);

        row = table(autocorrresults.scan_parameter_values(i), ...
            'VariableNames', {'scanParamVal'});
        summary_table = [summary_table; [row meanVals semVals]]; %#ok<AGROW>
    end

    %% --- Compute cumulants from full distribution of curves ---
    N_theta  = numel(autocorrresults.theta_values);
    g2_cums  = cell(N_params,1);

    for i = 1:N_params
        curves = autocorrresults.g2_curves{i}; % [N_reps × N_theta]
        cumul_mat = zeros(N_theta,4); % θ × first 4 cumulants
        for th = 1:N_theta
            g2vals = curves(:,th); % all repetitions at this theta
            cumul_mat(th,:) = Calculator.computeCumulants(g2vals,4);
        end
        g2_cums{i} = cumul_mat;
    end

    %% Package results
    results = struct();
    results.features_all    = feature_list;        % all repetitions pooled
    results.features_group  = per_group;           % per scan parameter
    results.meanCurve       = per_group_mean;      % mean-curve contrast & peaks
    results.summaries       = summary_table;       % group-level stats
    results.g2_cumulants    = g2_cums;             % cumulants from full distribution
    results.theta_values    = autocorrresults.theta_values;
end

function features = extractFeatures(thetas, g2curve)
    %% --- Step 1: Restrict theta to [0, π] ---
    mask           = (thetas >= 0 & thetas <= pi);
    thetasWithinPi = thetas(mask);
    g2WithinPi     = g2curve(mask);

    %% --- Step 2: DC removal ---
    g2cWithinPi    = g2WithinPi - mean(g2WithinPi);

    %% --- Step 3: Fourier projections (n=2,4,6) ---
    nlist = [2 4 6];
    An = zeros(size(nlist));
    for k = 1:numel(nlist)
        n = nlist(k);
        An(k) = abs(mean(g2cWithinPi .* exp(-1i*n*thetasWithinPi)));
    end

    %% --- Step 3b: Symmetry fractions ---
    totalAmp = sum(An);
    if totalAmp > 0
        S2 = An(1)/totalAmp;
        Q4 = An(2)/totalAmp;
        H6 = An(3)/totalAmp;
    else
        S2 = 0; Q4 = 0; H6 = 0;
    end

    %% --- Step 4: Contrast ---
    g2_restricted   = g2curve((thetas >= pi/18 & thetas <= 17*pi/18));
    
    if ~isempty(g2_restricted)
        ContrastVal = (max(g2_restricted) - min(g2_restricted)) / ...
                      (max(g2_restricted) + min(g2_restricted));
    else
        ContrastVal = NaN;
    end

    %% --- Step 8: Package as one row table ---
    features = table( ...
        An(1), An(2), An(3), ContrastVal, ...
        S2, Q4, H6, ...
        'VariableNames', {'A2','A4','A6','Contrast','S2','Q4','H6'});
end

function meanFeatures = extractMeanG2Features(thetas, g2_curves)
% Computes contrast and peak angles from the mean g2 curve.
% g2_curves: [N_reps × N_theta] array of individual g2(theta) curves
% thetas: vector of theta values in radians

    % --- Step 1: Compute mean curve ---
    meanG2    = mean(g2_curves, 1);

    % --- Step 2: Contrast calculation ---
    maskContrast    = (thetas >= pi/18 & thetas <= 17*pi/18);
    meanG2_contrast = meanG2(maskContrast);
    
    if ~isempty(meanG2_contrast)
        ContrastVal = (max(meanG2_contrast) - min(meanG2_contrast)) / ...
                      (max(meanG2_contrast) + min(meanG2_contrast));
    else
        ContrastVal = NaN;
    end

    % --- Step 3: Peak finding restricted to (pi/18, pi/2) ---
    maskPeaks    = (thetas >= pi/18 & thetas <= pi/2);
    meanG2_peaks = meanG2(maskPeaks);
    thetas_peaks = thetas(maskPeaks);

    % Dynamic prominence: 10% of range or minimum 0.05 absolute
    alpha       = 0.1;
    minAbsProm  = 0.005;
    dynamicProm = max(alpha * range(meanG2_peaks), minAbsProm);
    minDist     = max(round(numel(thetas_peaks)/12),1);  % at least 1

    [~, locs] = findpeaks(meanG2_peaks, ...
        'MinPeakProminence', dynamicProm, ...
        'MinPeakDistance', minDist);

    PeakAngles = thetas_peaks(locs);  % angles of detected peaks

    % --- Step 4: Fallback if no peaks detected ---
    if isempty(PeakAngles)
        % Use global maximum of restricted curve
        [~, idxMax] = max(meanG2_peaks);
        MaxPeakAngle = thetas_peaks(idxMax);
    else
        % Otherwise, take max of detected peaks
        [~, idxMax] = max(meanG2_peaks(locs));
        MaxPeakAngle = PeakAngles(idxMax);
    end

    % --- Step 5: Package as struct ---
    meanFeatures = struct();
    meanFeatures.MeanContrast    = ContrastVal;
    meanFeatures.MeanPeakAngles  = PeakAngles;
    meanFeatures.MaxPeakAngle    = MaxPeakAngle;

end
