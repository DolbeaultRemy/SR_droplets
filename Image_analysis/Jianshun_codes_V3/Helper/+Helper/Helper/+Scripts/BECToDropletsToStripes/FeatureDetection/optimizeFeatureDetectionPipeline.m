%% optimizePipeline.m
% Bayesian optimization of patch detection + shape extraction pipeline
% Joint optimization across all training images

% ------------------ USER INPUTS ------------------

% Use folder where this script is located
thisScriptPath        = mfilename('fullpath');
[thisScriptDir, ~, ~] = fileparts(thisScriptPath);
baseimageFolder       = fullfile(thisScriptDir, 'OptimizationImages');
basemaskFolder        = fullfile(thisScriptDir, 'OptimizationMasks');

% Load all images and ground truth masks
imgFiles  = dir(fullfile(baseimageFolder, '*.png'));
maskFiles = dir(fullfile(basemaskFolder, '*.png'));

nImages = min(numel(imgFiles), numel(maskFiles));
imgs    = cell(1, nImages);
masks   = cell(1, nImages);

for i = 1:nImages
    imgs{i}  = im2double(imread(fullfile(baseimageFolder, imgFiles(i).name)));
    masks{i} = imread(fullfile(basemaskFolder,  maskFiles(i).name)) > 0;
end

% ------------------ BASE PARAMETERS ------------------
params                           = struct();
params.backgroundDiskFraction    = 1/8;   % Fraction of image size used for morphological opening
params.boundingBoxPadding        = 12;    % Padding around detected cloud

% Initial (starting) values for optimizable parameters
params.dogGaussianSmallSigma     = 1.2412;
params.dogGaussianLargeSigma     = 3.9609;
params.adaptiveSensitivity       = 0.4146;
params.adaptiveNeighborhoodSize  = 17;
params.minPeakFraction           = 0.7610;
params.minimumPatchArea          = 80;
params.intensityThreshFraction   = 0.4879;
params.edgeSigma                 = 1.0244;
params.edgeThresholdLow          = 0.2384;
params.edgeThresholdHigh         = 0.6127;

% Fixed shape extraction parameters
params.shapeMinArea              = 20;
params.shapeCloseRadius          = 3;
params.shapeFillHoles            = false;

% Fixed imaging system metadata
params.pixelSize                 = 5.86e-6;   % meters/pixel
params.magnification             = 23.94;

% ------------------ OPTIMIZATION VARIABLES ------------------
optimVars = [
    optimizableVariable('dogGaussianSmallSigma',[1.0,1.5])
    optimizableVariable('dogGaussianLargeSigma',[3.5,4.0])
    optimizableVariable('adaptiveSensitivity',[0.3,0.6])
    optimizableVariable('adaptiveNeighborhoodSize',[11,19],'Type','integer')
    optimizableVariable('minPeakFraction',[0.7,0.85])
    optimizableVariable('minimumPatchArea',[80,150],'Type','integer')
    optimizableVariable('intensityThreshFraction',[0.35,0.55])
    optimizableVariable('edgeSigma',[0.8,1.2])
    optimizableVariable('edgeThresholdLow',[0.2,0.35])
    optimizableVariable('edgeThresholdHigh',[0.45,0.65])
    ];

% ------------------ OBJECTIVE FUNCTION ------------------
objFcn = @(x) mean(cellfun(@(im,mask) ...
                 Analyzer.evaluateFeatureDetectionPipeline(x, im, mask, params, true), ...
                 imgs, masks));

% ------------------ RUN BAYESIAN OPTIMIZATION ------------------
results = bayesopt(objFcn, ...
                   optimVars, ...
                   'MaxObjectiveEvaluations', 100, ...
                   'AcquisitionFunctionName','expected-improvement-plus', ...
                   'Verbose',1, ...
                   'PlotFcn',{@plotMinObjective,@plotAcquisitionFunction});

% ------------------ EXTRACT BEST PARAMETERS ------------------
bestX = results.XAtMinObjective;

params.dogGaussianSmallSigma     = bestX.dogGaussianSmallSigma;
params.dogGaussianLargeSigma     = bestX.dogGaussianLargeSigma;
params.adaptiveSensitivity       = bestX.adaptiveSensitivity;
params.adaptiveNeighborhoodSize  = bestX.adaptiveNeighborhoodSize;
params.minPeakFraction           = bestX.minPeakFraction;
params.minimumPatchArea          = bestX.minimumPatchArea;
params.intensityThreshFraction   = bestX.intensityThreshFraction;
params.edgeSigma                 = bestX.edgeSigma;
params.edgeThresholdLow          = bestX.edgeThresholdLow;
params.edgeThresholdHigh         = bestX.edgeThresholdHigh;

disp('Best parameters found (joint optimization across all images):');
disp(params);

%% --- Evaluate with best parameters, plot results on first image ---
fprintf('Evaluating with best parameters and plotting results...\n');

Analyzer.runInteractiveFeatureDetectorGUI(od_imgs, scan_parameter_values, file_list, options, params)
