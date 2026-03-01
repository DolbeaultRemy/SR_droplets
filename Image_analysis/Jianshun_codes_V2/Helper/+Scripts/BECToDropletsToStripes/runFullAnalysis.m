%% ===== BEC-Droplets-Stripes Settings =====

% Specify data location to run analysis on
dataSources = {
    struct('sequence', 'TwoDGas', ...
           'date', '2025/06/23', ...
           'runs', [300]) % specify run numbers as a string in "" or just as a numeric value
};

options = struct();

% File paths
options.baseDataFolder     = '//DyLabNAS/Data';
options.FullODImagesFolder = 'E:/Data - Experiment/FullODImages/202506';
options.measurementName    = 'DropletsToStripes';
scriptFullPath             = mfilename('fullpath');
options.saveDirectory      = fileparts(scriptFullPath);

% Camera / imaging settings
options.cam                = 4;             % 1 - ODT_1_Axis_Camera; 2 - ODT_2_Axis_Camera; 3 - Horizontal_Axis_Camera;, 4 - Vertical_Axis_Camera;
options.angle              = 0;             % angle by which image will be rotated
options.center             = [1410, 2030];
options.span               = [200, 200];
options.fraction           = [0.1, 0.1];
options.pixel_size         = 5.86e-6;       % in meters
options.magnification      = 23.94;
options.ImagingMode        = 'HighIntensity';
options.PulseDuration      = 5e-6;          % in s

% Fourier analysis settings
options.theta_min          = deg2rad(0);
options.theta_max          = deg2rad(180);
options.N_radial_bins      = 500;
options.Radial_Sigma       = 2;
options.Radial_WindowSize  = 5;              % odd number

options.k_min              = 1.2771;         % μm⁻¹
options.k_max              = 2.5541;         % μm⁻¹
options.N_angular_bins     = 180;
options.Angular_Threshold  = 75;
options.Angular_Sigma      = 2;
options.Angular_WindowSize = 5;
options.zoom_size          = 50;

% Flags
options.skipUnshuffling           = false;
options.skipNormalization         = false;

options.skipFringeRemoval         = true;
options.skipPreprocessing         = true;
options.skipMasking               = true;
options.skipIntensityThresholding = true;
options.skipBinarization          = true;

options.skipFullODImagesFolderUse = false;
options.skipSaveData              = false;
options.skipSaveFigures           = true;
options.skipSaveProcessedOD       = true;
options.skipLivePlot              = true;
options.showProgressBar           = true;

% Extras
options.font = 'Bahnschrift';
switch options.measurementName
    case 'BECToDroplets'
        options.scan_parameter        = 'rot_mag_field';
        options.flipSortOrder         = false;
        options.scanParameterUnits    = 'gauss';
        options.titleString           = 'BEC to Droplets';
    case 'BECToStripes'
        options.scan_parameter        = 'rot_mag_field';
        options.flipSortOrder         = false;
        options.scanParameterUnits    = 'gauss';
        options.titleString           = 'BEC to Stripes';
    case 'DropletsToStripes'
        options.scan_parameter        = 'ps_rot_mag_fin_pol_angle';
        options.flipSortOrder         = false;
        options.scanParameterUnits    = 'degrees';
        options.titleString           = 'Droplets to Stripes';
    case 'StripesToDroplets'
        options.scan_parameter        = 'ps_rot_mag_fin_pol_angle';
        options.flipSortOrder         = false;
        options.scanParameterUnits    = 'degrees';
        options.titleString           = 'Stripes to Droplets';
end

%% ===== Run Batch Analysis =====
results_all = Helper.batchAnalyze(dataSources, options);
