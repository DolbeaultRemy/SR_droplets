function results_all = batchAnalyze(dataSources, options)
%% batchAnalyze
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Brief description of the script functionality.
%
% Notes:
%   Optional notes, references.

    arguments
        dataSources (1,:) cell
        options struct
    end

    % Default base folder if not specified
    if ~isfield(options, 'baseDataFolder')
        options.baseDataFolder = '//DyLabNAS/Data';
    end

    % Default skip flag
    if ~isfield(options, 'skipFullODImagesFolderUse')
        options.skipFullODImagesFolderUse = false;
    end

    % Determine whether to use FullODImagesFolder or raw baseDataFolder
    useFullOD = false;
    if isfield(options, 'FullODImagesFolder') && isfolder(options.FullODImagesFolder)
        if ~isfolder(options.baseDataFolder)
            if ~options.skipFullODImagesFolderUse
                % Case 1a: raw folder missing, fallback to FullODImagesFolder
                useFullOD = true;
                fprintf('\n[INFO] Raw data folder not found but found full OD Images folder.\n');
                fprintf('\n[INFO] Using full OD Images folder: %s\n', options.FullODImagesFolder);
            else
                % Case 1b: raw folder missing, fallback to
                % FullODImagesFolder but user overrides
                error('Raw data folder not found, found full OD Images folder which cannot be used (set skipFullODImagesFolderUse=false to override). Aborting.\n');
            end
        elseif ~options.skipFullODImagesFolderUse
            % Case 2: both exist, prioritize FullODImagesFolder unless skipped
            useFullOD = true;
            fprintf('\n[INFO] Both raw data folder (%s) and full OD Images folder (%s) found.\n', ...
                options.baseDataFolder, options.FullODImagesFolder);
            fprintf('\n[INFO] Prioritizing full OD Images folder (set skipFullODImagesFolderUse=true to override).\n');
        else
            % Case 3: both exist, prioritize raw data folder because of overrride
            fprintf('\n[INFO] Both raw data folder (%s) and full OD Images folder (%s) found.\n', ...
                options.baseDataFolder, options.FullODImagesFolder);
            fprintf('\n[INFO] Prioritizing raw data folder (set skipFullODImagesFolderUse=false to override).\n');
        end
    elseif isfolder(options.baseDataFolder)
        % 🚨 Raw data folder exists, full OD images does not
        fprintf('\n[INFO] Full OD Images folder not found but found raw data folder.\n');
        fprintf('\n[INFO] Using raw data folder: %s\n', options.baseDataFolder);
    end

    % 🚨 Sanity check if neither folder exists
    if ~useFullOD && ~isfolder(options.baseDataFolder)
        warning('\n[ERROR] Neither raw data folder (%s) nor a valid full OD Images folder were found. Exiting.', ...
            options.baseDataFolder);
        results_all = {};
        return;
    end

    % ===== Estimate dataset memory and get per-run estimates =====
    [options.SAVE_TO_WORKSPACE, ~] = Helper.estimateDatasetMemory(dataSources, options);

    results_all = {};  % one element per folder

    % --- FULL OD MODE ---
    if useFullOD
        % --- List available FullODImages_* folders ---
        fullODFolders = dir(fullfile(options.FullODImagesFolder, 'FullODImages_*'));
        fullODFolders = fullODFolders([fullODFolders.isdir]);
    
        ds = dataSources{1}; % only one dataSources struct
    
        % Ensure sequences, dates, runs are cell arrays for uniform processing
        sequences = ds.sequence; if ischar(sequences), sequences = {sequences}; end
        dates     = ds.date;     if ischar(dates),     dates     = {dates};     end
        runs      = ds.runs;     if isnumeric(runs),   runs      = num2cell(runs); end
        if isstring(runs), runs = cellstr(runs); end
    
        % Loop over all combinations of sequence × date × run
        for seqIdx = 1:numel(sequences)
            for dateIdx = 1:numel(dates)
                for runIdx = 1:numel(runs)
                    targetSequence = sequences{seqIdx};
                    targetDate     = dates{dateIdx};
                    targetRun      = runs{runIdx};
    
                    matched = false;
                    for i = 1:numel(fullODFolders)
                        selectedPath = fullfile(fullODFolders(i).folder, fullODFolders(i).name);
    
                        % Load metadata for run info
                        metaFile = fullfile(selectedPath, 'metadata.mat');
                        if ~isfile(metaFile)
                            warning('No metadata.mat in %s, skipping.', selectedPath);
                            continue;
                        end
                        S = load(metaFile, 'metadata');
    
                        % Reconstruct sequence/date/run from stored folderPath
                        dataSourceMeta = makeDataSourceStruct(S.metadata.options.folderPath);
    
                        % Check for match: measurementName and run info
                        if isfield(S.metadata.options, 'measurementName') && ...
                           isfield(options, 'measurementName') && ...
                           strcmp(S.metadata.options.measurementName, options.measurementName) && ...
                           strcmp(dataSourceMeta{1}.sequence, targetSequence) && ...
                           strcmp(dataSourceMeta{1}.date, targetDate) && ...
                           dataSourceMeta{1}.runs == targetRun
    
                            fprintf('\n[INFO] Found matching full OD images subfolder: %s\n', fullODFolders(i).name);
                            options.selectedPath = selectedPath;
                            options.folderPath   = S.metadata.options.folderPath;
                            matched = true;
                            break;
                        end
                    end
    
                    if ~matched
                        warning('No matching full OD images subfolder found for sequence %s, date %s, run %s, measurementName %s.', ...
                               targetSequence, targetDate, targetRun, options.measurementName);
                        continue; % skip this run but continue to next combination
                    end
    
                    % ✅ Proceed to analysis for this combination
                    try
                        args = [fieldnames(options), struct2cell(options)]';
                        args = args(:)';
                        [analysisResults, scan_parameter_values, scan_reference_values] = Analyzer.performAnalysis(args{:});
    
                        result = struct();
                        result.sequence              = targetSequence;
                        result.date                  = targetDate;
                        result.run                   = targetRun;
                        result.path                  = options.folderPath;
                        result.options               = options;
                        result.results               = analysisResults;
                        result.scan_parameter_values = scan_parameter_values;
                        result.scan_reference_values = scan_reference_values;
    
                        % Save dataset as MAT
                        if ~isfield(options, 'skipSaveData') || ~options.skipSaveData
                            saveResultStruct(result, options.saveDirectory);
                        end
    
                        results_all{end+1,1} = result;
    
                    catch ME
                        warning("Error processing %s: %s", options.folderPath, ME.message);
                    end
                end
            end
        end
    
        return; % ✅ handled all in FullOD mode
    end

    % --- RAW MODE (default) ---
    ds = dataSources{1};  % single struct
    
    % Ensure sequences, dates, and runs are cell arrays
    sequences = ds.sequence; if ischar(sequences), sequences = {sequences}; end
    dates     = ds.date;     if ischar(dates),     dates     = {dates};     end
    runs      = ds.runs;     if isnumeric(runs),   runs      = num2cell(runs); end
    if isstring(runs), runs = cellstr(runs); end
    
    % Loop over all combinations of sequence × date × run
    for seqIdx = 1:numel(sequences)
        for dateIdx = 1:numel(dates)
            for runIdx = 1:numel(runs)
                targetSequence = sequences{seqIdx};
                targetDate     = dates{dateIdx};
                runItem        = runs{runIdx};
    
                % Convert runItem to string with leading zeros if numeric
                if isnumeric(runItem)
                    runID = sprintf('%04d', runItem);
                elseif isstring(runItem)
                    runID = runItem;
                elseif ischar(runItem)
                    runID = string(runItem);
                elseif iscell(runItem)
                    runID = string(runItem{1});
                else
                    error('Unsupported type for run entry: %s', class(runItem));
                end
    
                % Determine base folder
                if isfield(ds, 'baseFolder') && ~isempty(ds.baseFolder)
                    baseFolder = fullfile(ds.baseFolder, targetSequence, targetDate);
                else
                    baseFolder = fullfile(options.baseDataFolder, targetSequence, targetDate);
                end
                
                % Build folder path
                folderPath         = fullfile(baseFolder, runID);
                options.folderPath = folderPath;
    
                try
                    % Convert struct -> name-value args
                    args = [fieldnames(options), struct2cell(options)]';
                    args = args(:)';
    
                    % Perform analysis
                    [analysisResults, scan_parameter_values] = Analyzer.performAnalysis(args{:});
    
                    % Store results
                    result = struct();
                    result.sequence              = targetSequence;
                    result.date                  = targetDate;
                    result.run                   = runID;
                    result.path                  = folderPath;
                    result.options               = options;
                    result.results               = analysisResults;
                    result.scan_parameter_values = scan_parameter_values;
    
                    % Save dataset
                    if ~isfield(options, 'skipSaveData') || ~options.skipSaveData
                        saveResultStruct(result, options.saveDirectory);
                    end
    
                    % Append to output
                    results_all{end+1,1} = result;
    
                catch ME
                    warning("Error processing %s/%s/%s: %s", ...
                            targetSequence, targetDate, runID, ME.message);
                end
            end
        end
    end
end

%% --- Local helper functions ---
function saveResultStruct(result, saveDirectory)
    % Define results folder
    resultsFolder = fullfile(saveDirectory, "Results", "SavedData");
    if ~exist(resultsFolder, 'dir')
        mkdir(resultsFolder);
    end

    % Path to index file
    indexFile = fullfile(resultsFolder, "datasetsIndex.mat");

    % Load or initialize index
    if isfile(indexFile)
        S = load(indexFile, "nextIdx");
        nextIdx = S.nextIdx;
    else
        nextIdx = 1;
    end

    % Variable name and file path
    varName = sprintf('Dataset_%d', nextIdx);
    savePath = fullfile(resultsFolder, varName + ".mat");

    % Save dataset as struct inside MAT file
    S.(varName) = result;
    save(savePath, '-struct', 'S');

    % Update index
    nextIdx = nextIdx + 1;
    save(indexFile, "nextIdx");
end

function dataSource = makeDataSourceStruct(folderPath)
    % Split by file separators (handles / or \)
    parts = regexp(folderPath, '[\\/]', 'split');

    % Remove empty parts caused by leading slashes
    parts = parts(~cellfun('isempty', parts));

    % Extract sequence, date, and run number
    % Now the indices are correct:
    % parts = {'DyLabNAS', 'Data', 'StructuralPhaseTransition', '2025', '08', '13', '0062'}
    sequence = parts{3};       % "StructuralPhaseTransition"
    year     = parts{4};       % "2025"
    month    = parts{5};       % "08"
    day      = parts{6};       % "13"
    runStr   = parts{7};       % "0062"

    % Build date string
    dateStr = sprintf('%s/%s/%s', year, month, day);

    % Convert run string to number
    runNum = str2double(runStr);

    % Construct struct inside a cell array
    dataSource = {
        struct('sequence', sequence, ...
               'date', dateStr, ...
               'runs', runNum)
    };
end

