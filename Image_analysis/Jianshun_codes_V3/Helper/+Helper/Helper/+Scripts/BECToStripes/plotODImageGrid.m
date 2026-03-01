%% JULY ===== BEC-Stripes Settings =====

% Specify data location to run analysis on
dataSources = {
    struct('sequence', 'StructuralPhaseTransition', ...
           'date', '2025/07/26', ...
           'runs', [8]) % specify run numbers as a string in "" or just as a numeric value
};

options = struct();

% File paths
options.baseDataFolder     = '//DyLabNAS/Data';
options.FullODImagesFolder = 'E:/Data - Experiment/FullODImages/202507';
options.measurementName    = 'BECToStripes';
scriptFullPath             = mfilename('fullpath');
options.saveDirectory      = fileparts(scriptFullPath);

% Camera / imaging settings
options.cam                = 4;             % 1 - ODT_1_Axis_Camera; 2 - ODT_2_Axis_Camera; 3 - Horizontal_Axis_Camera;, 4 - Vertical_Axis_Camera;
options.angle              = 0;             % angle by which image will be rotated
options.center             = [1435, 2035];
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
        options.flipSortOrder         = true;
        options.scanParameterUnits    = 'degrees';
        options.titleString           = 'Droplets to Stripes';
    case 'StripesToDroplets'
        options.scan_parameter        = 'ps_rot_mag_fin_pol_angles';
        options.flipSortOrder         = false;
        options.scanParameterUnits    = 'degrees';
        options.titleString           = 'Stripes to Droplets';
end

% ===== Collect Images and Launch Viewer =====

[options.selectedPath, options.folderPath]                 = Helper.selectDataSourcePath(dataSources, options);

[od_imgs, scan_parameter_values, scan_reference_values, ~] = Helper.collectODImages(options);

% === Plot raw images ===

% Select every 10th repetition
selected_reps              = [1, 10, 20, 30, 40, 50];
nRepsPlot                  = numel(selected_reps);

% === Select subset of scan parameter values (x-values) to plot ===
% Example: skip some B-fields
xValsToPlot                = [2.5, 2.45, 2.40, 2.35, 2.20, 2.00];  % <-- keep only these

[~, idxKeep]               = ismembertol(xValsToPlot, scan_reference_values, 1e-6);

missing                    = xValsToPlot(idxKeep==0);
if ~isempty(missing)
    warning('The following values were not found in scan_reference_values: %s', num2str(missing));
end
idxKeep                    = idxKeep(idxKeep>0);

scan_reference_values_plot = scan_reference_values(idxKeep);
nParamsPlot                = numel(scan_reference_values_plot);

% === Collect ALL repetitions for each kept parameter ===
od_imgs_grouped            = cell(nParamsPlot, 1);

for k = 1:nParamsPlot
    paramVal = scan_reference_values_plot(k);

    % Find all repetitions for this parameter value
    repIdx = find(ismembertol(scan_parameter_values, paramVal, 1e-6));

    % Store images + parameter values
    od_imgs_grouped{k}     = od_imgs(repIdx);
end

figure(101); clf;
set(gcf, 'Position', [100 100 950 750]);

% Rows = repetitions, Cols = selected scan parameter values
t = tiledlayout(nRepsPlot, nParamsPlot, 'TileSpacing', 'compact', 'Padding', 'compact');

font = 'Bahnschrift';
allAxes = gobjects(nRepsPlot, nParamsPlot);

for r = 1:nRepsPlot
    for j = 1:nParamsPlot
        ax = nexttile((r-1)*nParamsPlot + j);
        allAxes(r,j) = ax;

        % Extract OD image for this scan parameter and repetition
        od = od_imgs_grouped{j}{r};

        imagesc(od, 'Parent', ax);
        set(ax, 'YDir', 'normal');
        axis(ax, 'image');
        ax.DataAspectRatioMode = 'manual';
        ax.PlotBoxAspectRatioMode = 'auto';

        ax.XTick = [];
        ax.YTick = [];

        colormap(ax, Colormaps.inferno());
    end
end

% Shared colorbar
cb = colorbar('Location', 'eastoutside');
ylabel(cb, 'Optical Density');
cb.Layout.Tile = 'east';
cb.FontName = font;
cb.FontSize = 18;
cb.Label.FontSize = 20;
cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.HorizontalAlignment = 'center';
cb.Direction = 'normal';

% Label B-field values along bottom
for j = 1:nParamsPlot
    ax = allAxes(end, j);
    ax.XTick = size(od,2)/2;
    ax.XTickLabel = sprintf('%.2f G', scan_reference_values_plot(j));
    ax.XTickLabelRotation = 45;
    ax.FontName = font;
    ax.FontSize = 16;
end

% Label repetitions along left
for r = 1:nRepsPlot
    ax = allAxes(r, 1);
    ax.YTick = size(od,1)/2;
    ax.YTickLabel = sprintf('Run %d', selected_reps(r));
    ax.FontName = font;
    ax.FontSize = 16;
end

% === Add figure title ===
title(t, sprintf('%s: \\alpha = %.2f^\\circ', options.titleString, 180-145), ...
    'FontSize', 24, 'FontName', font, 'FontWeight', 'bold');


%% AUGUST ===== BEC-Stripes Settings =====

% Specify data location to run analysis on
dataSources = {
    struct('sequence', 'StructuralPhaseTransition', ...
           'date', '2025/08/15', ...
           'runs', [3]) % specify run numbers as a string in "" or just as a numeric value
};

options = struct();

% File paths
options.baseDataFolder     = '//DyLabNAS/Data';
options.FullODImagesFolder = 'E:/Data - Experiment/FullODImages/202508';
options.measurementName    = 'BECToStripes';
scriptFullPath             = mfilename('fullpath');
options.saveDirectory      = fileparts(scriptFullPath);

% Camera / imaging settings
options.cam                = 4;             % 1 - ODT_1_Axis_Camera; 2 - ODT_2_Axis_Camera; 3 - Horizontal_Axis_Camera;, 4 - Vertical_Axis_Camera;
options.angle              = 0;             % angle by which image will be rotated
options.center             = [1420, 2050];
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
        options.flipSortOrder         = true;
        options.scanParameterUnits    = 'degrees';
        options.titleString           = 'Droplets to Stripes';
    case 'StripesToDroplets'
        options.scan_parameter        = 'ps_rot_mag_fin_pol_angles';
        options.flipSortOrder         = false;
        options.scanParameterUnits    = 'degrees';
        options.titleString           = 'Stripes to Droplets';
end

% ===== Collect Images and Launch Viewer =====

[options.selectedPath, options.folderPath]                 = Helper.selectDataSourcePath(dataSources, options);

[od_imgs, scan_parameter_values, scan_reference_values, ~] = Helper.collectODImages(options);

% === Plot raw images ===
n_total                    = numel(scan_parameter_values);
nParams                    = numel(scan_reference_values);  
nReps                      = n_total/nParams;

% Select every 10th repetition
selected_reps              = [1, 10, 20, 30, 40, 50];
nRepsPlot                  = numel(selected_reps);

% === Select subset of scan parameter values (x-values) to plot ===
% Example: skip some B-fields
xValsToPlot                = [2.45, 2.40, 2.35, 2.30, 2.20, 2.10, 2.00];  % <-- keep only these

[~, idxKeep]               = ismembertol(xValsToPlot, scan_reference_values, 1e-6);

missing                    = xValsToPlot(idxKeep==0);
if ~isempty(missing)
    warning('The following values were not found in scan_reference_values: %s', num2str(missing));
end
idxKeep                    = idxKeep(idxKeep>0);

scan_reference_values_plot = scan_reference_values(idxKeep);
nParamsPlot                = numel(scan_reference_values_plot);

% === Collect ALL repetitions for each kept parameter ===
od_imgs_grouped            = cell(nParamsPlot, 1);

for k = 1:nParamsPlot
    paramVal = scan_reference_values_plot(k);

    % Find all repetitions for this parameter value
    repIdx = find(ismembertol(scan_parameter_values, paramVal, 1e-6));

    % Store images + parameter values
    od_imgs_grouped{k}     = od_imgs(repIdx);
end

figure(102); clf;
set(gcf, 'Position', [100 100 950 750]);

% Rows = repetitions, Cols = selected scan parameter values
t = tiledlayout(nRepsPlot, nParamsPlot, 'TileSpacing', 'compact', 'Padding', 'compact');

font = 'Bahnschrift';
allAxes = gobjects(nRepsPlot, nParamsPlot);

for r = 1:nRepsPlot
    for j = 1:nParamsPlot
        ax = nexttile((r-1)*nParamsPlot + j);
        allAxes(r,j) = ax;

        % Extract OD image for this scan parameter and repetition
        od = od_imgs_grouped{j}{r};

        imagesc(od, 'Parent', ax);
        set(ax, 'YDir', 'normal');
        axis(ax, 'image');
        ax.DataAspectRatioMode = 'manual';
        ax.PlotBoxAspectRatioMode = 'auto';

        ax.XTick = [];
        ax.YTick = [];

        colormap(ax, Colormaps.inferno());
    end
end

% Shared colorbar
cb = colorbar('Location', 'eastoutside');
ylabel(cb, 'Optical Density');
cb.Layout.Tile = 'east';
cb.FontName = font;
cb.FontSize = 18;
cb.Label.FontSize = 20;
cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.HorizontalAlignment = 'center';
cb.Direction = 'normal';

% Label B-field values along bottom
for j = 1:nParamsPlot
    ax = allAxes(end, j);
    ax.XTick = size(od,2)/2;
    ax.XTickLabel = sprintf('%.2f G', scan_reference_values_plot(j));
    ax.XTickLabelRotation = 45;
    ax.FontName = font;
    ax.FontSize = 16;
end

% Label repetitions along left
for r = 1:nRepsPlot
    ax = allAxes(r, 1);
    ax.YTick = size(od,1)/2;
    ax.YTickLabel = sprintf('Run %d', selected_reps(r));
    ax.FontName = font;
    ax.FontSize = 16;
end

% === Add figure title ===
title(t, sprintf('%s: \\alpha = %.2f^\\circ', options.titleString, 180-145), ...
    'FontSize', 24, 'FontName', font, 'FontWeight', 'bold');
