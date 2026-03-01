%% ------------------ 1. Compare radial spectral contrast across datasets ------------------
datasetIdxList = 1:6;

% --- Base directory selection ---
useLocalBaseDir = false;  % <-- set true to use script location, false to use manual path

if useLocalBaseDir
    % Use folder where this script is located
    thisScriptPath        = mfilename('fullpath');
    [thisScriptDir, ~, ~] = fileparts(thisScriptPath);
    baseDir               = fullfile(thisScriptDir, 'Results');
else
    % Use manually specified folder
    baseDir               = 'E:\Results - Experiment\202508\BECToDropletsToStripes\';
end

% Prepare storage
scanValsCell   = cell(1, numel(datasetIdxList));
meanValsCell   = cell(1, numel(datasetIdxList));
stderrValsCell = cell(1, numel(datasetIdxList));
labelsCell     = {'2.25 G', '2.20 G', '2.15 G', '2.10 G', '2.05 G', '2.00 G'};

% --- Load options from first dataset to set plotting defaults ---
datasetIdx  = datasetIdxList(1);
datasetName = sprintf('Dataset_%d', datasetIdx);
dataFile    = fullfile(baseDir, 'SavedData', [datasetName '.mat']);
data        = load(dataFile);
datasetStruct = data.(datasetName);
options = datasetStruct.options;
options.font            = 'Bahnschrift';
options.skipSaveFigures = false;

% --- Ensure combined figure folder exists ---
combinedSaveDir = fullfile(baseDir, 'SavedFigures', 'Combined');
if ~exist(combinedSaveDir, 'dir')
    mkdir(combinedSaveDir);
end

% --- Loop over datasets to extract mean ± SE ---
for i = 1:numel(datasetIdxList)
    datasetIdx  = datasetIdxList(i);
    datasetName = sprintf('Dataset_%d', datasetIdx);

    % Build paths for this dataset
    dataFile   = fullfile(baseDir, 'SavedData', [datasetName '.mat']);
    
    % Load dataset
    data          = load(dataFile);
    datasetStruct = data.(datasetName);

    % Extract values
    scanVals = datasetStruct.scan_parameter_values;
    dataVals = datasetStruct.results.spectral_analysis_results.radial_spectral_contrast;

    % --- Compute mean and standard error ---
    [unique_vals, ~, idx] = unique(scanVals);
    mean_vals   = zeros(size(unique_vals));
    stderr_vals = zeros(size(unique_vals));
    for k = 1:length(unique_vals)
        if iscell(dataVals)
            group = dataVals{idx == k};
        else
            group = dataVals(idx == k);
        end
        if iscell(group)
            groupVals = [group{:}];
        else
            groupVals = group;
        end
        mean_vals(k)   = mean(groupVals);
        stderr_vals(k) = std(groupVals) / sqrt(length(groupVals));
    end

    % Store results
    scanValsCell{i}   = unique_vals;
    meanValsCell{i}   = mean_vals;
    stderrValsCell{i} = stderr_vals;
end

% --- Call compare function ---
Plotter.compareMultipleDatasets(scanValsCell, meanValsCell, stderrValsCell, ...
    'FigNum', 16, ...
    'FontName', options.font, ...
    'Labels', labelsCell, ...
    'Title', 'Radial Spectral Contrast Across Datasets', ...
    'XLabel', 'B (G)', ...
    'YLabel', 'Radial Spectral Contrast', ...
    'SaveDirectory', combinedSaveDir, ...
    'SaveFileName', 'RadialSpectralContrast_Combined.fig', ...
    'SkipSaveFigures', false);

%% ------------------ 2. Compare max g2 across datasets ------------------

% Prepare storage
maxG2ValsCell   = cell(1, numel(datasetIdxList));
maxG2MeanCell   = cell(1, numel(datasetIdxList));
maxG2StderrCell = cell(1, numel(datasetIdxList));

for i = 1:numel(datasetIdxList)
    datasetIdx  = datasetIdxList(i);
    datasetName = sprintf('Dataset_%d', datasetIdx);

    % Load dataset
    dataFile      = fullfile(baseDir, 'SavedData', [datasetName '.mat']);
    data          = load(dataFile);
    datasetStruct = data.(datasetName);

    % Extract values
    scanVals = datasetStruct.scan_parameter_values;
    dataVals = datasetStruct.results.custom_g_results.max_g2_all_per_scan_parameter_value;

    % --- Compute mean and standard error ---
    [unique_vals, ~, idx] = unique(scanVals);
    mean_vals   = zeros(size(unique_vals));
    stderr_vals = zeros(size(unique_vals));
    for k = 1:length(unique_vals)
        if iscell(dataVals)
            group = dataVals{idx == k};
        else
            group = dataVals(idx == k);
        end
        if iscell(group)
            groupVals = [group{:}];
        else
            groupVals = group;
        end
        mean_vals(k)   = mean(groupVals);
        stderr_vals(k) = std(groupVals) / sqrt(length(groupVals));
    end

    % Store results
    maxG2ValsCell{i}   = unique_vals;
    maxG2MeanCell{i}   = mean_vals;
    maxG2StderrCell{i} = stderr_vals;
end

% --- Call compare function ---
Plotter.compareMultipleDatasets(maxG2ValsCell, maxG2MeanCell, maxG2StderrCell, ...
    'FigNum', 17, ...
    'FontName', options.font, ...
    'Labels', labelsCell, ...
    'Title', 'Peak Offset Angular Correlation Across Datasets', ...
    'XLabel', 'B (G)', ...
    'YLabel', '$\mathrm{max}[g^{(2)}_{[50,70]}(\delta\theta)]$', ...
    'SaveDirectory', combinedSaveDir, ...
    'SaveFileName', 'MaxG2_Combined.fig', ...
    'SkipSaveFigures', false);