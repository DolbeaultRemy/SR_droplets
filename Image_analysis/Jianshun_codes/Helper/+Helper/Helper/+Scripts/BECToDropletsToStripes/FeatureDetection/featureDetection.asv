%% ===== BEC-Droplets-Stripes Settings =====

% Specify data location to run analysis on
dataSources = {
    struct('sequence', 'StructuralPhaseTransition', ...
           'date', '2025/11/03', ...
           'runs', "0006") % specify run numbers as a string in "" or just as a numeric value
};

options = struct();

% File paths
options.baseDataFolder     = '//DyLabNAS/Data';
options.FullODImagesFolder = 'E:/Data - Experiment/FullODImages/202506';
options.measurementName    = 'DropletsToStripes';%'DropletsToStripes';
scriptFullPath             = mfilename('fullpath');
options.saveDirectory      = fileparts(scriptFullPath);

% Camera / imaging settings
options.cam                = 4;             % 1 - ODT_1_Axis_Camera; 2 - ODT_2_Axis_Camera; 3 - Horizontal_Axis_Camera;, 4 - Vertical_Axis_Camera;
options.angle              = 0;             % angle by which image will be rotated
options.center             = [1410, 2030];
options.span               = [200, 200];
options.fraction           = [0.1, 0.1];
options.pixel_size         = 5.86e-6;       % in meters
options.magnification      = 24.6;
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
        options.flipSortOrder         = true;
        options.scanParameterUnits    = 'gauss';
        options.titleString           = 'BEC to Droplets';
    case 'BECToStripes'
        options.scan_parameter        = 'rot_mag_field';
        options.flipSortOrder         = true;
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

%% ===== Collect Images and Launch Viewer =====

[options.selectedPath, options.folderPath]                         = Helper.selectDataSourcePath(dataSources, options);

[od_imgs, scan_parameter_values, scan_reference_values, file_list] = Helper.collectODImages(options);

%%
save_folder = 'C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2025\11\03\0006'; % Windows路径示例
% save_folder = '/home/username/your/folder/path/'; % Linux/Mac路径示例

% 确保文件夹存在，如果不存在则创建
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

num_images = length(od_imgs);

for i = 1:num_images
    % 生成文件名，使用4位数字编号（前导零）
    filename = sprintf('Image_%04d.mat', i);
    
    % 完整的文件路径
    full_path = fullfile(save_folder, filename);
    
    % 获取当前图片
    OD = od_imgs{i};
    
    % 保存为.mat文件
    save(full_path, 'OD', 'scan_parameter_values', 'scan_reference_values', 'scan_parameter_names');
    
    % 可选：显示进度
    if mod(i, 50) == 0
        fprintf('已保存 %d/%d 个图片\n', i, num_images);
    end
end

%% Minimal pipeline: preprocess + detect patches + extract shapes

% ---------- user-tunable parameters ----------
params.backgroundDiskFraction    = 0.1250;          % Fraction of image size used for morphological opening
params.boundingBoxPadding        = 12;              % Padding around detected cloud
params.dogGaussianSmallSigma     = 1.1498;          % Sigma for small Gaussian in Difference-of-Gaussians
params.dogGaussianLargeSigma     = 3.9232;          % Sigma for large Gaussian in Difference-of-Gaussians
params.adaptiveSensitivity       = 0.4552;          % Adaptive threshold sensitivity Higher → more pixels marked foreground (lower threshold).
params.adaptiveNeighborhoodSize  = 13;              % Window size for adaptive threshold Defines the local window size over which threshold is computed. Larger → smoother masks, less sensitive to noise.
params.minPeakFraction           = 0.7190;          % Fraction of max DoG response to reject weak patches
params.minimumPatchArea          = 81;              % Minimum area (pixels) for detected patches to be kept

params.shapeMinArea              = 20;              % Minimum shape area
params.shapeCloseRadius          = 3;               % Morphological closing radius (fills holes)
params.shapeFillHoles            = false;           % Ensures shapes are solid regions
params.intensityThreshFraction   = 0.4499;          % Fraction of max intensity to keep
params.edgeSigma                 = 1.1749;          % Gaussian smoothing for Canny
params.edgeThresholdLow          = 0.3383;          % Low Canny threshold fraction
params.edgeThresholdHigh         = 0.6412;          % High Canny threshold fraction

params.pixelSize                 = 5.86e-6;         % Physical pixel size (meters/pixel)
params.magnification             = 23.94;           % Magnification factor of imaging system
% --------------------------------------------

Analyzer.runInteractiveFeatureDetectorGUI(od_imgs, scan_parameter_values, file_list, options, params)


