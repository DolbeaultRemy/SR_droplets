function runAnalyse(dataSources, folderPath, ROIcenter, measurementName)
% RUNANALYSE Main function to process experimental data.
%   dataSources: cell array of structs with fields 'sequence', 'date', 'runs'
%   folderPath: string, path where processed data will be saved
%   ROIcenter: [x y] center coordinates for cropping (optional)
%   measurementName: string, e.g., 'BECToDroplets' (optional)

    %% ----- Set default parameters and options -----
    if nargin < 3 || isempty(ROIcenter)
        ROIcenter = []; % will be set later based on measurement type if needed
    end
    if nargin < 4 || isempty(measurementName)
        measurementName = 'Default';
    end

    options = configureOptions(measurementName, ROIcenter, folderPath);

    %% ----- Feature detection parameters -----
    detectionParams = setDetectionParameters();

    % DEBUGGING FLAG – set to true to save data before detectStructure
    debugDetect = true;   % ← change to false when not debugging

    %% ----- Collect images and scan parameters -----
    [selectedPath, options.folderPath] = Helper.selectDataSourcePath(dataSources, options);
    [odImages, scanParamValues, scanRefValues, fileList] = Helper.collectODImages(options);
    scanParamNames = options.scanParameter;

    %% ----- Save raw OD images as .mat files -----
    ensureDirectoryExists(folderPath);
    numImages = length(odImages);
    for i = 1:numImages
        filename = sprintf('Image_%04d.mat', i);
        fullPath = fullfile(folderPath, filename);
        OD = odImages{i};
        save(fullPath, 'OD', 'scanParamValues', 'scanRefValues', 'scanParamNames');
        if mod(i, 50) == 0
            fprintf('Saved %d/%d images\n', i, numImages);
        end
    end

    %% ----- Feature detection parameters -----
    detectionParams = setDetectionParameters();

    %% ----- Process each image for structure detection and bottleneck analysis -----
    matFiles = dir(fullfile(folderPath, '*.mat'));
    for i = 1:length(matFiles)
        load(fullfile(folderPath, matFiles(i).name));
        img = double(OD);

        % Optionally crop the image if ROIcenter is provided
        % if ~isempty(options.center)
        %     img = Helper.cropODImage(img, options.center, options.span);
        % end

        % ---------- DEBUG: save inputs to detectStructure ----------
        if debugDetect
            debugFolder = fullfile(folderPath, 'debug_detectStructure');
            if ~exist(debugFolder, 'dir')
                mkdir(debugFolder);
            end
            debugFileName = fullfile(debugFolder, ['debug_' matFiles(i).name]);
            % Save the image, parameters, and some context
            save(debugFileName, 'img', 'detectionParams', 'options', 'matFiles', 'i');
            fprintf('Debug data saved to %s\n', debugFileName);
        end
        % ------------------------------------------------------------

        % Detect structures (patches)
        [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart, binaryMask, CC] = ...
            detectStructure(img, detectionParams, folderPath, matFiles(i).name);

        % Perform bottleneck analysis and splitting
        [bottleneckResults, bwSplit, ccSplit, splitInfo] = ...
            comprehensiveBottleneckAnalysisWithSplitting(binaryMask, CC, 0.40, 10);

        % Visualize results
        visualizeComprehensiveBottleneckAnalysisWithSplitting(binaryMask, ccSplit, ...
            bottleneckResults, imgCropped, folderPath, matFiles(i).name);

        % Save analysis data
        saveAnalysisData(folderPath, matFiles(i).name, ccSplit, imgCropped, ...
            scanParamValues, scanRefValues, scanParamNames);
    end
end

%% ----- Helper functions -----
function options = configureOptions(measurementName, ROIcenter, folderPath)
% Configure options struct based on measurement type.
    options = struct();
    options.baseDataFolder     = '//DyLabNAS/Data';
    options.FullODImagesFolder = 'C:\Users\Jianshun Gao\Documents\DyData';
    options.measurementName    = measurementName;
    scriptFullPath             = mfilename('fullpath');
    options.saveDirectory      = fileparts(scriptFullPath);

    % Camera settings
    options.cam                = 4;
    options.angle              = 0;
    options.center             = ROIcenter;
    options.span               = [300, 300];
    options.fraction           = [0.1, 0.1];
    options.pixel_size         = 5.86e-6;
    options.magnification      = 24.6;
    options.ImagingMode        = 'HighIntensity';
    options.PulseDuration      = 5e-6;

    % Fourier analysis settings
    options.theta_min          = deg2rad(0);
    options.theta_max          = deg2rad(180);
    options.N_radial_bins      = 500;
    options.Radial_Sigma       = 2;
    options.Radial_WindowSize  = 5;
    options.k_min              = 1.2771;
    options.k_max              = 2.5541;
    options.N_angular_bins     = 180;
    options.Angular_Threshold  = 75;
    options.Angular_Sigma      = 2;
    options.Angular_WindowSize = 5;
    options.zoom_size          = 50;

    % Processing flags
    options.skipUnshuffling           = false;
    options.skipNormalization         = false;
    options.skipFringeRemoval         = true;
    options.skipPreprocessing         = true;
    options.skipMasking               = true;
    options.skipIntensityThresholding = true;
    options.skipBinarization          = true;
    options.skipFullODImagesFolderUse = true;
    options.skipSaveData              = true;
    options.skipSaveFigures           = true;
    options.skipSaveProcessedOD       = true;
    options.skipLivePlot              = true;
    options.showProgressBar           = true;

    % Measurement-specific parameters
    options.font = 'Bahnschrift';
    switch measurementName
        case 'BECToDroplets'
            options.scanParameter     = 'rot_mag_field';
            options.flipSortOrder      = true;
            options.scanParameterUnits = 'gauss';
            options.titleString        = 'BEC to Droplets';
        case 'DropletsToBEC'
            options.scanParameter     = 'ps_rot_mag_field';
            options.flipSortOrder      = true;
            options.scanParameterUnits = 'gauss';
            options.titleString        = 'Droplets to BEC';
        case 'BECToStripes'
            options.scanParameter     = 'rot_mag_field';
            options.flipSortOrder      = true;
            options.scanParameterUnits = 'gauss';
            options.titleString        = 'BEC to Stripes';
        case 'DropletsToStripes'
            options.scanParameter     = 'rot_mag_field';
            options.flipSortOrder      = false;
            options.scanParameterUnits = 'degrees';
            options.titleString        = 'Droplets to Stripes';
        case 'StripesToDroplets'
            options.scanParameter     = 'ps_rot_mag_fin_pol_angle';
            options.flipSortOrder      = false;
            options.scanParameterUnits = 'degrees';
            options.titleString        = 'Stripes to Droplets';
        otherwise
            options.scanParameter     = 'rot_mag_field';
            options.flipSortOrder      = false;
            options.scanParameterUnits = 'a.u.';
            options.titleString        = measurementName;
    end
end

function params = setDetectionParameters()
% Set parameters for structure detection.
    params.backgroundDiskFraction    = 0.1250;
    params.boundingBoxPadding        = 35;
    params.dogGaussianSmallSigma     = 0.5;
    params.dogGaussianLargeSigma     = 4;
    params.adaptiveSensitivity       = 0.3;
    params.adaptiveNeighborhoodSize  = 13;
    params.minPeakFraction           = 0.2;
    params.minimumPatchArea          = 20;
    params.shapeMinArea              = 20;
    params.shapeCloseRadius          = 3;
    params.shapeFillHoles            = false;
    params.intensityThreshFraction   = 0.4499;
    params.edgeSigma                 = 1.1749;
    params.edgeThresholdLow          = 0.3383;
    params.edgeThresholdHigh         = 0.6412;
    params.pixelSize                 = 5.86e-6;
    params.magnification             = 23.94;
end

function ensureDirectoryExists(folderPath)
% Create directory if it does not exist.
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
end

function saveAnalysisData(folderPath, fileName, ccSplit, imgCropped, ...
    scanParamValues, scanRefValues, scanParamNames)
% Save the analysis results to a .mat file.
    savePath = fullfile(folderPath, 'AnalysisPlot', 'Data');
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    save(fullfile(savePath, fileName), 'ccSplit', 'imgCropped', ...
        'scanParamValues', 'scanRefValues', 'scanParamNames');
end