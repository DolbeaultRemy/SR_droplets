%% ===== BEC-Droplets-Stripes Settings =====

% Specify data location to run analysis on
dataSources = {
    struct('sequence', 'StructuralPhaseTransition', ...
           'date', '2025/08/20', ...
           'runs', [1]) % specify run numbers as a string in "" or just as a numeric value
};

options = struct();

% File paths
options.baseDataFolder     = '//DyLabNAS/Data';
options.FullODImagesFolder = 'E:/Data - Experiment/FullODImages/202508';
options.measurementName    = 'DropletsToStripes';
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
options.skipLivePlot              = false;
options.showProgressBar           = true;

% Extras
options.font = 'Bahnschrift';
switch options.measurementName
    case 'DropletsToStripes'
        options.ignore_scan_parameter = {'ps_rot_mag_field'}; % Parameter to IGNORE from these - rot_mag_field, ps_rot_mag_field, ps_rot_mag_fin_pol_angle
        options.flipSortOrder         = true;
        options.scanParameterUnits    = {'gauss','degrees'};
        options.titleString           = 'Droplets to Stripes';
    case 'StripesToDroplets'
        options.ignore_scan_parameter = {'ps_rot_mag_field'}; % Parameter to IGNORE from these - rot_mag_field, ps_rot_mag_field, ps_rot_mag_fin_pol_angle
        options.flipSortOrder         = false;
        options.scanParameterUnits    = {'gauss','degrees'};
        options.titleString           = 'Stripes to Droplets';
end

%% ===== Collect Images and Launch Viewer =====

[options.selectedPath, options.folderPath]                 = Helper.selectDataSourcePath(dataSources, options);

[od_imgs, scan_parameter_values, scan_reference_values, ~] = Helper.collectODImages(options);

%% === Plot raw images ===

% === User selection of BFields and alpha values ===
BFieldsToPlot = [2.45, 2.4, 2.35, 2.34, 2.32, 2.3, 2.28];  % [2.45, 2.44, 2.43, 2.42, 2.4, 2.39, 2.38, 2.37, 2.36, 2.35, 2.34, 2.32, 2.3, 2.28, 2.25, 2.2, 2.15, 2.1, 2.05, 2.0, 1.95, 1.90, 1.85, 1.8] G
alphasToPlot  = [0, 5, 10, 20, 30, 40];  % [0, 5, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40] degrees

% --- Convert scan_reference_values (cell array of 1x2) into Nx2 array ---
scanRefArray     = cell2mat(scan_reference_values(:));  % Nx2
scanParamArray   = cell2mat(scan_parameter_values(:));  % Nx2

BFieldVals    = scanRefArray(:,1);
alphaVals     = scanRefArray(:,2);

% --- Select subset of scan_reference_values to plot ---
keepIdx       = ismembertol(BFieldVals, BFieldsToPlot, 1e-6) & ...
                ismembertol(alphaVals, alphasToPlot, 1e-6);

scan_reference_values_plot = scan_reference_values(keepIdx);
nParamsPlot                = numel(scan_reference_values_plot);

% --- Warn about missing values ---
existingB     = BFieldsToPlot(ismembertol(BFieldsToPlot, unique(BFieldVals), 1e-6));
missingB      = setdiff(BFieldsToPlot, existingB);

existingA     = alphasToPlot(ismembertol(alphasToPlot, unique(alphaVals), 1e-6));
missingA      = setdiff(alphasToPlot, existingA);

if ~isempty(missingB)
    warning('The following BFields were not found: %s', num2str(missingB));
end
if ~isempty(missingA)
    warning('The following alpha values were not found: %s', num2str(missingA));
end

BFieldsPlot   = existingB;
alphasPlot    = existingA;

% --- Collect all repetitions for each kept parameter ---
od_imgs_grouped = cell(nParamsPlot, 1);

for k = 1:nParamsPlot
    paramVal = scan_reference_values_plot{k};  % 1x2 vector [BField, alpha]

    % Find all repetitions matching this parameter pair (BField, alpha)
    repIdx = find(ismembertol(scanParamArray, paramVal, 1e-6, 'ByRows', true));

    % Store images + parameter values
    od_imgs_grouped{k} = od_imgs(repIdx);
end

% === Plot OD images as a phase diagram ===
% Unique values for axes
BFieldsPlot = sort(unique(cellfun(@(c) c(1), scan_reference_values_plot)), 'descend'); % y-axis
alphasPlot  = sort(unique(cellfun(@(c) c(2), scan_reference_values_plot)));            % x-axis

nB = numel(BFieldsPlot);
nA = numel(alphasPlot);

repToPlot = 1;  % which repetition to show

% Create figure
figure(300); clf;
set(gcf, 'Position', [100 100 950 750]);
t = tiledlayout(nB, nA, 'TileSpacing', 'compact', 'Padding', 'compact');

font = 'Bahnschrift';
allAxes = gobjects(nB, nA);

% Loop over rows (BFields) and columns (alphas)
for r = 1:nB
    B = BFieldsPlot(r);
    for c = 1:nA
        alpha = alphasPlot(c);

        ax = nexttile((r-1)*nA + c);
        allAxes(r,c) = ax;

        % Find the index k corresponding to this (BField, alpha) pair
        k = find(cellfun(@(v) all(abs(v - [B alpha]) < 1e-6), scan_reference_values_plot), 1);

        if ~isempty(k)
            img = od_imgs_grouped{k}{repToPlot};
            imagesc(ax, img);
        else
            % If missing, leave blank
            axis(ax, 'off');
            continue
        end

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

% Label alpha values along bottom
for c = 1:nA
    ax = allAxes(end, c);
    ax.XTick = size(img,2)/2;
    ax.XTickLabel = sprintf('%.0f°', alphasPlot(c));
    ax.XTickLabelRotation = 45;
    ax.FontName = font;
    ax.FontSize = 16;
end

% Label BFields along left
for r = 1:nB
    ax = allAxes(r, 1);
    ax.YTick = size(img,1)/2;
    ax.YTickLabel = sprintf('%.2f G', BFieldsPlot(r));
    ax.FontName = font;
    ax.FontSize = 16;
end

% Figure title
title(t, sprintf('%s', options.titleString), ...
    'FontSize', 24, 'FontName', font, 'FontWeight', 'bold');