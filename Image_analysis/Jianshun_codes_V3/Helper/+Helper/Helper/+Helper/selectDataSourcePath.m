function [selectedPath, folderPath] = selectDataSourcePath(dataSources, options)
%% selectDataSourcePath
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Helper function to select a run/folder for analysis.
%
% Inputs:
%   dataSources - cell array of structs with fields: sequence, date, runs
%   options     - struct, may contain: baseDataFolder, FullODImagesFolder, skipFullODImagesFolderUse, saveDirectory
%
% Outputs:
%   selectedPath - actual folder to use (raw or FullOD)
%   folderPath   - constructed path from dataSources (always raw folder)
%
% Notes:
%   Optional notes, references.
    
    allPaths = {};  % initialize candidate paths
    lookup   = struct('rawPath', {}, 'fullODPath', {});
    
    % --- General path to full od image folders ---
    if isfield(options, 'FullODImagesFolder')
        full_od_image_parent_folder = options.FullODImagesFolder;
    elseif isfield(options, 'saveDirectory')
        full_od_image_parent_folder = options.saveDirectory;
    else
        full_od_image_parent_folder = '';
    end

    % --- Gather candidate raw data paths ---
    for i = 1:numel(dataSources)
        ds = dataSources{i};

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
                        runID = char(runItem);
                    elseif ischar(runItem)
                        runID = runItem;
                    elseif iscell(runItem)
                        runID = char(runItem{1});
                    else
                        error('Unsupported type for run entry: %s', class(runItem));
                    end

                    % Build raw data path
                    dateParts = strsplit(targetDate,'/');
                    rawPath = fullfile(options.baseDataFolder, targetSequence, dateParts{:}, runID);
                    if isfolder(rawPath)
                        allPaths{end+1} = rawPath;
                    else
                        fprintf('\n[INFO] Raw data folder does not exist: %s\n', rawPath);
                    end

                    % Build matching FullOD path
                    expectedName = sprintf('FullODImages_%s_%s_Run%s', ...
                                           targetSequence, strrep(targetDate,'/','-'), runID);
                    if isfolder(full_od_image_parent_folder) && ~isempty(full_od_image_parent_folder)
                        fullODPath = fullfile(full_od_image_parent_folder, expectedName);
                    else
                        fullODPath = '';
                    end

                    % Store lookup mapping
                    lookup(end+1).rawPath   = rawPath; %#ok<AGROW>
                    lookup(end).fullODPath  = fullODPath;
                end
            end
        end
    end

    % --- Determine whether FullODImagesFolder should be used ---
    useFullOD = false;
    if ~isempty(allPaths)
        if isfolder(full_od_image_parent_folder) && ~isempty(full_od_image_parent_folder)
            if ~isfield(options,'skipFullODImagesFolderUse') || ~options.skipFullODImagesFolderUse
                fprintf('\n[INFO] Both raw data folder (%s) and full OD Images folder (%s) found.\n', ...
                    options.baseDataFolder, full_od_image_parent_folder);
                fprintf('\n[INFO] Prioritizing full OD Images folder (set skipFullODImagesFolderUse=true to override).\n');
                useFullOD = true;
            else
                fprintf('\n[INFO] Both raw data folder (%s) and full OD Images folder (%s) found.\n', ...
                    options.baseDataFolder, full_od_image_parent_folder);
                fprintf('\n[INFO] Prioritizing raw data folder (set skipFullODImagesFolderUse=false to override).\n');
            end
        else
            fprintf('\n[INFO] Using raw data folder(s) since full OD images not found or not specified.\n');
        end
    elseif isfolder(full_od_image_parent_folder) && ~isempty(full_od_image_parent_folder)
        if ~isfield(options,'skipFullODImagesFolderUse') || ~options.skipFullODImagesFolderUse
            useFullOD = true;
            fprintf('\n[INFO] Raw data folder(s) not found but found full OD Images folder which will be used.\n');
            allPaths{end+1} = full_od_image_parent_folder;
        else
            error('Raw data folder(s) not found, found full OD Images folder which cannot be used (set skipFullODImagesFolderUse=false to override). Aborting.\n');
        end
    end

    if isempty(allPaths)
        error('No valid paths for data found. Aborting.');
    end

    % --- If using full OD images folder, match the subfolder corresponding to sequence/date/run ---
    if useFullOD
        matchedPaths = {};
        for k = 1:numel(lookup)
            if isfolder(lookup(k).fullODPath)
                matchedPaths{end+1} = lookup(k).fullODPath;
            end
        end
        if ~isempty(matchedPaths)
            allPaths = matchedPaths;
        else
            error('No valid FullODImages_* subfolders found under %s. Aborting.', full_od_image_parent_folder);
        end
    end

    % --- Determine if user selection is needed ---
    if numel(allPaths) > 1
        % Build compact display names for selection
        listStrings = cell(size(allPaths));
        for idx = 1:numel(allPaths)
            % Match to lookup table
            label = '';
            for k = 1:numel(lookup)
                if strcmp(allPaths{idx}, lookup(k).rawPath) || strcmp(allPaths{idx}, lookup(k).fullODPath)
                    % Extract from rawPath parts
                    rawParts = strsplit(lookup(k).rawPath, filesep);
                    runID    = rawParts{end};
                    seq      = rawParts{end-4};
                    year     = rawParts{end-3};
                    month    = rawParts{end-2};
                    day      = rawParts{end-1};
                    dateFormatted = sprintf('%s/%s/%s', day, month, year); % dd/mm/yyyy
                    label = sprintf('%s | %s | Run:%s', seq, dateFormatted, runID);
                    break;
                end
            end
            if isempty(label)
                [~, label] = fileparts(allPaths{idx}); % fallback
            end
            listStrings{idx} = label;
        end
    
        set(0,'DefaultUicontrolFontSize',10);
        [selectedIndex, tf] = listdlg('PromptString','Select a run to analyze:', ...
                                      'SelectionMode','single', ...
                                      'ListString', listStrings, ...
                                      'ListSize',[500, 300]);
        if ~tf
            error('No path selected. Aborting.');
        end
        selectedPath = allPaths{selectedIndex};
    else
        selectedPath = allPaths{1};
    end

    % --- Construct corresponding raw folderPath ---
    folderPath = '';
    for k = 1:numel(lookup)
        if strcmp(selectedPath, lookup(k).rawPath) || strcmp(selectedPath, lookup(k).fullODPath)
            folderPath = lookup(k).rawPath;
            break;
        end
    end

    if isempty(folderPath)
        error('Could not map selected path back to a raw data folder path.');
    end

    fprintf('\n[INFO] Selected path: %s\n', selectedPath);

end
