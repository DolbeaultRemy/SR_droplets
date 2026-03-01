function [od_imgs, scan_parameter_values, scan_reference_values, file_list] = collectODImages(options)
%% collectODImages
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Applies cropping, background subtraction, and optional fringe removal, optional unshuffling on OD image dataset.
%   Automatically reuses in-memory full dataset if available otherwise, reads and processes raw HDF5 data.
%
% Inputs:
%   options - structure containing processing options:
%       .folderPath           : path to raw HDF5 files
%       .saveDirectory        : path to save cache (if needed)
%       .cam, .angle          : camera selection and rotation angle
%       .ImagingMode, .PulseDuration : imaging parameters
%       .scan_parameter       : name of scan parameter
%       .center, .span        : cropping settings
%       .fraction             : background subtraction fraction
%       .removeFringes        : logical flag for fringe removal
%       .skipUnshuffling      : logical flag to skip unshuffling
%       .scan_reference_values: reference values for unshuffling
%
% Outputs:
%        od_imgs              : cell array of processed OD images
%        scan_parameter_values: vector of scan parameter values
%        file_list            : cell array of file names
%
% Notes:
%   Optional notes, references.

    % --- Early exit if processed OD images already exist AND options match ---
    reuseExistingODImages = evalin('base', ...
        'exist(''od_imgs'',''var'') && exist(''scan_parameter_values'',''var'') && exist(''file_list'',''var'') && exist(''prior_options'',''var'')');

    if reuseExistingODImages
        prior_options = evalin('base','prior_options');
        critical_fields = {'folderPath','cam','angle','ImagingMode','PulseDuration',...
                           'center','span','fraction','removeFringes','skipUnshuffling','scan_reference_values'};

        if ~haveOptionsChanged(options, prior_options, critical_fields)
            fprintf('\n[INFO] Reusing processed OD images, scan parameters, and file list from memory.\n');
            od_imgs               = evalin('base','od_imgs');
            scan_parameter_values = evalin('base','scan_parameter_values');
            scan_reference_values = evalin('base','scan_reference_values');
            file_list             = evalin('base','file_list');

            % --- Save OD images to disk if requested ---
            if ~options.skipSaveProcessedOD
                saveProcessedOD(od_imgs, options);
            end

            return; % ✅ bypass cropping/background subtraction
        else
            fprintf('\n[INFO] Processed-data-related options changed. Reprocessing full OD image dataset...\n');
        end
    end

    % --- General path to full od image folders ---
    if isfield(options, 'FullODImagesFolder')
        full_od_image_parent_folder = options.FullODImagesFolder;
    elseif isfield(options, 'saveDirectory')
        full_od_image_parent_folder = options.saveDirectory;
    else
        full_od_image_parent_folder = '';
    end

    % --- Specific sequence, data and run ---
    dataSource       = makeDataSourceStruct(options.folderPath);

    % --- Check if workspace full dataset exists ---
    fullDataExists = evalin('base', 'exist(''full_od_imgs'', ''var'')') && ...
                     evalin('base', 'exist(''full_bkg_imgs'', ''var'')') && ...
                     evalin('base', 'exist(''raw_scan_parameter_values'', ''var'')') && ...
                     evalin('base', 'exist(''raw_file_list'', ''var'')') && ...
                     evalin('base', 'exist(''prior_options'',''var'')');

    if ~isfield(options,'SAVE_TO_WORKSPACE')
        [options.SAVE_TO_WORKSPACE, ~] = Helper.estimateDatasetMemory(dataSource, options);
    end
    
    full_od_image_subfolder               = [];

    % --- Prepare full_od_imgs, full_bkg_imgs, scan values, file list ---
    if fullDataExists
        % --- Case 1: Already in workspace ---
        fprintf('\n[INFO] Reusing full OD image dataset and scan parameters from memory.\n');
        full_od_imgs              = evalin('base','full_od_imgs');
        full_bkg_imgs             = evalin('base','full_bkg_imgs');
        raw_scan_parameter_values = evalin('base','raw_scan_parameter_values');
        raw_file_list             = evalin('base','raw_file_list');
        nFiles                    = size(full_od_imgs,3);
        fprintf('\n[INFO] Cropping and subtracting background from images...\n');
    
    else
        matched = false;
        useFullODFolders   = ~isfield(options,'skipFullODImagesFolderUse') || ~options.skipFullODImagesFolderUse;
    
        % --- If user provided a selected path and it is a folder, disambiguate its meaning ---
        if isfield(options,'selectedPath') && isfolder(options.selectedPath)
            selPath = options.selectedPath;
    
            % --- Determine if selectedPath looks like a FullODImages folder ---
            selIsFullOD = false;
            if isfile(fullfile(selPath,'metadata.mat'))
                selIsFullOD = true;
            else
                % exact match to any discovered fullodimage_folders entry
                if isfolder(full_od_image_parent_folder) && ~isempty(full_od_image_parent_folder)
                    for r = 1:numel(full_od_image_parent_folder)
                        if strcmp(full_od_image_parent_folder, selPath)
                            selIsFullOD = true;
                            break;
                        end
                    end
                end
            end
    
            if selIsFullOD && useFullODFolders
                % --- selectedPath is explicitly the full OD images folder ---
                full_od_image_subfolder = selPath;
                matched = true;
                fprintf('\n[INFO] Using selected full OD images subfolder: %s\n', full_od_image_subfolder);
    
            else
                % --- selectedPath appears to be a raw-data folder (not a full-OD folder) ---
                fprintf('\n[INFO] Selected path appears to be raw data: %s\n', selPath);
    
                % If user forces recompute -> recompute from selected raw path
                if isfield(options,'skipFullODImagesFolderUse') && options.skipFullODImagesFolderUse
                    [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, raw_scan_parameter_names, scan_reference_values] = recomputeODImages(options, selPath);
                    if ~options.SAVE_TO_WORKSPACE
                        % --- Determine parent folder for FullODImages ---
                        if isfield(options, 'FullODImagesFolder') && ...
                           isfolder(options.FullODImagesFolder) && ...
                           ~isempty(options.FullODImagesFolder)
                            parentFolder = options.FullODImagesFolder;
                        else
                            parentFolder = options.saveDirectory;
                        end
                        dataSource                = makeDataSourceStruct(options.folderPath);
                        full_od_image_subfolder   = createFullODImagesFolderPath(parentFolder, dataSource);
                        useFullODFolders          = true;
                    else
                        nFiles  = numel(raw_file_list);
                    end
                    matched = true;    
                else
                    % Try to find an existing full-OD folder whose metadata references this raw path
                    found = false;
                    if isfolder(full_od_image_parent_folder) && ~isempty(full_od_image_parent_folder) && useFullODFolders
                        for r = 1:numel(full_od_image_parent_folder)
                            metaPath = fullfile(full_od_image_parent_folder(r).folder, full_od_image_parent_folder(r).name, 'metadata.mat');
                            if ~isfile(metaPath), continue; end
                            S = load(metaPath,'metadata');
                            % Compare data source from metadata to the selected raw path
                            mdDataSource = makeDataSourceStruct(S.metadata.options.folderPath);
                            selDataSource = makeDataSourceStruct(selPath);
                            if isfield(S.metadata.options,'measurementName') && isfield(options,'measurementName') && ...
                               strcmp(S.metadata.options.measurementName, options.measurementName) && ...
                               isequal(mdDataSource, selDataSource)
                                full_od_image_subfolder = fullfile(full_od_image_parent_folder(r).folder, full_od_image_parent_folder(r).name);
                                matched = true;
                                found = true;
                                fprintf('\n[INFO] Found matching full OD images subfolder for selected raw path: %s\n', full_od_image_parent_folder(r).name);
                                break;
                            end
                        end
                    end
    
                    % If no matching full-OD folder found, recompute from the selected raw path
                    if ~found
                        fprintf('\n[INFO] Forcing recompute from raw data as no matching full OD images subfolder found (Use skipFullODImagesFolderUse=true to skip directly to computing from raw data).\n');
                        [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, raw_scan_parameter_names, scan_reference_values] = recomputeODImages(options, selPath);
                        if ~options.SAVE_TO_WORKSPACE
                            % --- Determine parent folder for FullODImages ---
                            if isfield(options, 'FullODImagesFolder') && ...
                               isfolder(options.FullODImagesFolder) && ...
                               ~isempty(options.FullODImagesFolder)
                                parentFolder = options.FullODImagesFolder;
                            elseif isfield(options, 'saveDirectory') && isfolder(options.saveDirectory)
                                parentFolder = options.saveDirectory;
                            end
                            dataSource                = makeDataSourceStruct(options.folderPath);
                            full_od_image_subfolder   = createFullODImagesFolderPath(parentFolder, dataSource);
                            useFullODFolders          = true;
                        else
                            nFiles  = numel(raw_file_list);
                        end
                        matched = true;
                    end
                end
            end
    
        else
            % --- No selected path: either force recompute or search among fullodimage_folders ---
            if isfield(options,'skipFullODImagesFolderUse') && options.skipFullODImagesFolderUse
                [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, raw_scan_parameter_names, scan_reference_values] = recomputeODImages(options, options.baseDataFolder);
                if ~options.SAVE_TO_WORKSPACE
                    % --- Determine parent folder for FullODImages ---
                    if isfield(options, 'FullODImagesFolder') && ...
                       isfolder(options.FullODImagesFolder) && ...
                       ~isempty(options.FullODImagesFolder)
                        parentFolder = options.FullODImagesFolder;
                    elseif isfield(options, 'saveDirectory') && isfolder(options.saveDirectory)
                        parentFolder = options.saveDirectory;
                    end
                    dataSource                = makeDataSourceStruct(options.folderPath);
                    full_od_image_subfolder   = createFullODImagesFolderPath(parentFolder, dataSource);
                    useFullODFolders          = true;
                else
                    nFiles  = numel(raw_file_list);
                end
                matched = true;
            else
                % Search for existing matching full-OD folder based on options.folderPath
                if isfolder(full_od_image_parent_folder) && ~isempty(full_od_image_parent_folder) && useFullODFolders
                    for r = 1:numel(full_od_image_parent_folder)
                        metaPath = fullfile(full_od_image_parent_folder(r).folder, full_od_image_parent_folder(r).name,'metadata.mat');
                        if ~isfile(metaPath), continue; end
                        S = load(metaPath,'metadata');
    
                        mdDataSource      = makeDataSourceStruct(S.metadata.options.folderPath);
                        currentDataSource = makeDataSourceStruct(options.folderPath);
    
                        if isfield(S.metadata.options,'measurementName') && isfield(options,'measurementName') && ...
                           strcmp(S.metadata.options.measurementName, options.measurementName) && ...
                           isequal(mdDataSource, currentDataSource)
                            full_od_image_subfolder = fullfile(full_od_image_parent_folder(r).folder, full_od_image_parent_folder(r).name);
                            matched = true;
                            fprintf('\n[INFO] Found matching full OD images subfolder: %s\n', full_od_image_parent_folder(r).name);
                            break;
                        end
                    end
                end
            end
        end
    
        % --- If still not matched, recompute from raw (fallback) ---
        if ~matched
            if ~isfolder(full_od_image_parent_folder) && useFullODFolders
                fprintf('\n[INFO] No full OD images found in workspace or on disk. Will recompute from raw data.\n');
            elseif isfolder(full_od_image_parent_folder) && useFullODFolders
                fprintf('\n[INFO] No matching full OD images subfolder found. Will recompute from raw data.\n');
            end
            [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, raw_scan_parameter_names, scan_reference_values] = recomputeODImages(options, options.baseDataFolder);
            if ~options.SAVE_TO_WORKSPACE
                % --- Determine parent folder for FullODImages ---
                if isfield(options, 'FullODImagesFolder') && ...
                   isfolder(options.FullODImagesFolder) && ...
                   ~isempty(options.FullODImagesFolder)
                    parentFolder = options.FullODImagesFolder;
                elseif isfield(options, 'saveDirectory') && isfolder(options.saveDirectory)
                    parentFolder = options.saveDirectory;
                end
                dataSource                = makeDataSourceStruct(options.folderPath);
                full_od_image_subfolder   = createFullODImagesFolderPath(parentFolder, dataSource);
                useFullODFolders          = true;
            else
                nFiles  = numel(raw_file_list);
            end
        end
    
        % --- If a folder was determined, load its contents (listing) ---
        if ~isempty(full_od_image_subfolder) && useFullODFolders
            if isfolder(full_od_image_subfolder)
                [mat_files, raw_scan_parameter_names, raw_scan_parameter_values, raw_file_list, nFiles] = prepareFromOnDiskData(full_od_image_subfolder);
                fprintf('\n[INFO] Cropping and subtracting background from images in full OD images folder on disk...\n');
            end
        end
    end
    
    % --- Unified cropping & background subtraction ---
    absimages = zeros(options.span(1)+1, options.span(2)+1, nFiles, 'single');
    refimages = zeros(options.span(1)+1, options.span(2)+1, nFiles, 'single');

    showPB = isfield(options,'showProgressBar') && options.showProgressBar;
    if showPB
        pb = Helper.ProgressBar();
        pb.run('Progress: ');
    end

    for k = 1:nFiles
        if fullDataExists || ~isempty(full_od_image_subfolder)
            if fullDataExists
                od_mat  = full_od_imgs(:,:,k);
                bkg_mat = full_bkg_imgs(:,:,k);
            else
                data                         = load(fullfile(full_od_image_subfolder, mat_files(k).name));
                od_mat                       = data.OD;
                bkg_mat                      = data.BKG;
                
                % --- Handle parameter names ---
                raw_scan_parameter_names{k} = data.Scan_Param;
    
                % --- Handle parameter values ---
                if numel(data.Scan_Val) == 1
                    % First time through: initialize numeric storage if empty
                    if isempty(raw_scan_parameter_values)
                        raw_scan_parameter_values = zeros(1,nFiles);
                    elseif iscell(raw_scan_parameter_values)
                        error('Mixed single-parameter and multi-parameter scans detected.');
                    end
                    raw_scan_parameter_values(k) = data.Scan_Val;
    
                else
                    % First time through: initialize cell storage if empty
                    if isempty(raw_scan_parameter_values)
                        raw_scan_parameter_values = cell(nFiles,1);
                    elseif isnumeric(raw_scan_parameter_values)
                        error('Mixed single-parameter and multi-parameter scans detected.');
                    end
                    raw_scan_parameter_values{k} = data.Scan_Val(:).';
                end

                raw_file_list(k)             = data.File;
            end
        else
            od_mat  = full_od_imgs(:,:,k);
            bkg_mat = full_bkg_imgs(:,:,k);
        end

        % --- Crop and subtract background ---
        if any(isnan(od_mat(:)))
            absimages(:,:,k) = nan(options.span(1)+1, options.span(2)+1, 'single');
        else
            cropped_od = Helper.cropODImage(od_mat, options.center, options.span);
            absimages(:,:,k) = Helper.subtractBackgroundOffset(cropped_od, options.fraction)';
        end

        if any(isnan(bkg_mat(:)))
            refimages(:,:,k) = nan(options.span(1)+1, options.span(2)+1, 'single');
        else
            cropped_bkg = Helper.cropODImage(bkg_mat, options.center, options.span);
            refimages(:,:,k) = Helper.subtractBackgroundOffset(cropped_bkg, options.fraction)';
        end

        if showPB
            progressPercent = round(k/nFiles*100);
            pb.run(progressPercent);
        end
    end

    if showPB
        pb.run(' Done!');
    end

    % --- Optional fringe removal ---
    if isfield(options, 'skipFringeRemoval') && ~options.skipFringeRemoval
        fprintf('\n[INFO] Applying fringe removal to processed images...\n');
        optrefimages             = Helper.removeFringesInImage(absimages, refimages);
        absimages_fringe_removed = absimages - optrefimages;
        od_imgs                  = arrayfun(@(i) absimages_fringe_removed(:,:,i), 1:size(absimages,3), 'UniformOutput', false);
        fprintf('\n[INFO] Fringe removal completed.\n');
    else
        od_imgs                  = arrayfun(@(i) absimages(:,:,i), 1:size(absimages,3), 'UniformOutput', false);
    end
    
    % --- Optional unshuffling based on scan reference values ---
    if isfield(options, 'skipUnshuffling') && ~options.skipUnshuffling
        fprintf('\n[INFO] Reordering images...\n');
    
        % --- Determine scan reference values ---
        if ~isfield(options, 'scan_reference_values') || isempty(options.scan_reference_values)
            if isnumeric(raw_scan_parameter_values) && isvector(raw_scan_parameter_values)
                % --- Single parameter case (numeric vector) ---
                scan_reference_values = unique(raw_scan_parameter_values(:), 'stable');
                n_total               = numel(raw_scan_parameter_values);
                n_values              = numel(scan_reference_values);
            elseif iscell(raw_scan_parameter_values)
                % --- Multi-parameter case (cell array of row vectors) ---
                params                = cell2mat(raw_scan_parameter_values); % convert to numeric matrix
                scan_reference_values = unique(params, 'rows', 'stable');
                n_total               = numel(raw_scan_parameter_values);
                n_values              = size(scan_reference_values, 1);
            else
                error('Unsupported format for raw scan parameter values.');
            end
        else
            % --- Use reference values from options ---
            scan_reference_values = options.scan_reference_values;
        
            if isnumeric(raw_scan_parameter_values) && isvector(raw_scan_parameter_values)
                n_total  = numel(raw_scan_parameter_values);
                n_values = numel(scan_reference_values);
            else
                n_total  = size(raw_scan_parameter_values, 1);
                n_values = size(scan_reference_values, 1);
            end
        end
        
        % --- Determine sort order ---
        if isfield(options, 'flipSortOrder') && options.flipSortOrder
            sort_order = 'descend';
        else
            sort_order = 'ascend'; % default
        end
        
        % --- Sort reference values ---
        if isnumeric(scan_reference_values) && isvector(scan_reference_values)
            scan_reference_values = sort(scan_reference_values, sort_order);
        else
            scan_reference_values = sortrows(scan_reference_values, sort_order);
        end
        
        % --- Convert multi-parameter scan_reference_values to cell array ---
        if ~isnumeric(scan_reference_values) || size(scan_reference_values,2) > 1
            scan_reference_values = mat2cell(scan_reference_values, ...
                                             ones(size(scan_reference_values,1),1), ...
                                             size(scan_reference_values,2));
        end

        % --- Reorder images according to scan reference values ---
        n_reps   = n_total / n_values;
        tol      = 1e-6;  % tolerance for floating-point comparisons
        idx_mat  = nan(n_reps, n_values);
        
        for j = 1:n_values
            if isnumeric(raw_scan_parameter_values)
                % Single parameter
                ref_val = scan_reference_values(j);  % scalar
                idx_all = find(abs(raw_scan_parameter_values - ref_val) < tol);
            elseif iscell(raw_scan_parameter_values)
                % Multi-parameter
                ref_val = scan_reference_values{j};    % row vector
                diffs   = cellfun(@(x) all(abs(x - ref_val) < tol), raw_scan_parameter_values);
                idx_all = find(diffs);
            end
        
            if numel(idx_all) ~= n_reps
                error('Reference value(s) %s occurs %d times; expected %d', ...
                      mat2str(ref_val,6), numel(idx_all), n_reps);
            end
        
            idx_mat(:, j) = idx_all(:);
        end
        
        % --- Reordered indices ---
        ordered_idx = reshape(idx_mat.', 1, []);
        
        % --- Apply reorder ---
        od_imgs = od_imgs(ordered_idx);
        
        if isnumeric(raw_scan_parameter_values)
            scan_parameter_values = raw_scan_parameter_values(ordered_idx).';
        elseif iscell(raw_scan_parameter_values)
            scan_parameter_values = raw_scan_parameter_values(ordered_idx);
        end
        
        file_list = raw_file_list(ordered_idx);
        
        fprintf('\n[INFO] Image reordering completed.\n');
    
    else
        scan_parameter_values = raw_scan_parameter_values;
        file_list             = raw_file_list;
    end
    
    % --- Determine scan parameter(s) ---
    if ~isfield(options,'scan_parameter') || isempty(options.scan_parameter)
        
        % Flatten all names into a single cell array
        all_names = {};
        for k = 1:numel(raw_scan_parameter_names)
            x = raw_scan_parameter_names{k};
            if iscell(x)
                all_names = [all_names, x(:).'];  % flatten row
            else
                all_names{end+1} = x;             % single char array
            end
        end
        
        % Find unique names (stable order)
        unique_names = unique(all_names, 'stable');
        
        % Decide single vs multiple parameter output
        if numel(unique_names) == 1
            scan_parameter_names   = unique_names{1}; % single char array
        else
            scan_parameter_names   = unique_names;        % cell array of char arrays
        end
    else
        scan_parameter_names   = options.scan_parameter;
    end
    
    % --- Save processed dataset and options to workspace ---
    assignin('base', 'od_imgs', od_imgs);
    assignin('base', 'scan_parameter_names', scan_parameter_names);
    assignin('base', 'scan_parameter_values', scan_parameter_values);
    assignin('base', 'scan_reference_values', scan_reference_values);
    assignin('base', 'file_list', file_list);
    assignin('base', 'prior_options', options);
    
    % --- Save OD images to disk if requested ---
    if ~options.skipSaveProcessedOD
        saveProcessedOD(od_imgs, options);
    end
    
    fprintf('\n[INFO] OD image dataset ready for further analysis.\n');
end

%% --- Local helper functions ---

function changed = haveOptionsChanged(options, prior_options, critical_fields)
    changed = false;
    for f = critical_fields
        fname = f{1};
        if isfield(options, fname) && isfield(prior_options, fname)
            if ~isequal(options.(fname), prior_options.(fname))
                changed = true; return
            end
        elseif xor(isfield(options, fname), isfield(prior_options, fname))
            changed = true; return
        end
    end
end

function saveProcessedOD(od_imgs, options)
% Saves a cell array of processed OD images (with metadata) to .mat files
% od_imgs is expected to be a cell array of structs with fields: OD, Scan, File
%
% Inputs:
%   od_imgs - cell array of structs: OD (2D array), Scan (scalar), File (string)
%   options - struct containing folderPath and either FullODImagesFolder or saveDirectory

    nImgs = numel(od_imgs);
    if nImgs == 0
        error('[ERROR] No images found.');
    end

    % --- Determine image size ---
    [ny, nx] = size(od_imgs{1}.OD);

    % --- Create unique folder name ---
    dataSource = makeDataSourceStruct(options.folderPath);
    runID = sprintf('%s_%s_Run%04d', ...
                    dataSource{1}.sequence, ...
                    strrep(dataSource{1}.date,'/','-'), ...
                    dataSource{1}.runs);

    % --- Determine parent folder ---
    if isfield(options, 'FullODImagesFolder') && isfolder(options.FullODImagesFolder) && ~isempty(options.FullODImagesFolder)
        parentFolder = options.FullODImagesFolder;
    elseif isfield(options, 'saveDirectory') && isfolder(options.saveDirectory)
        parentFolder = options.saveDirectory;
    else
        parentFolder = pwd;  % fallback to current folder
    end

    % --- Create ProcessedODImages folder ---
    processedFolder = fullfile(parentFolder, ['ProcessedODImages_' runID]);
    if ~exist(processedFolder, 'dir')
        mkdir(processedFolder);
    end
    fprintf('\n[INFO] Saving processed OD images in folder: %s\n', processedFolder);

    % --- Check if files already exist ---
    filesExist = all(arrayfun(@(k) isfile(fullfile(processedFolder, sprintf('Image_%04d.mat', k))), 1:nImgs)) ...
                 && isfile(fullfile(processedFolder,'metadata.mat'));

    if filesExist
        fprintf('\n[INFO] Processed OD .mat files already exist in %s. Skipping save.\n', processedFolder);
        return;
    end

    % --- Save metadata ---
    metadata.options   = options;
    metadata.timestamp = datetime;
    metadata.runID     = runID;
    metadata.imageSize = [ny, nx];
    metadata.fileList  = string(cellfun(@(c) c.File, od_imgs, 'UniformOutput', false));
    save(fullfile(processedFolder,'metadata.mat'), 'metadata', '-v7.3');

    % --- Save each image as a struct with OD, Scan, File ---
    for k = 1:nImgs
        imgStruct = od_imgs{k};
        if ~isstruct(imgStruct) || ~all(isfield(imgStruct, {'OD','Scan','File'}))
            error('od_imgs{%d} must be a struct with fields OD, Scan, File.', k);
        end

        OD   = single(imgStruct.OD);
        Scan = single(imgStruct.Scan);
        File = string(imgStruct.File);

        matFilePath = fullfile(processedFolder, sprintf('Image_%04d.mat', k));
        save(matFilePath, 'OD','Scan','File','-v7.3');
    end

    fprintf('\n[INFO] Processed OD .mat files and metadata saved successfully.\n');
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

function [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, scan_parameter_names, scan_reference_values] = recomputeODImages(options, selectedPath)
    % recomputeODImages: central recompute routine
    % optional argument selectedPath (if provided) will be used as options.folderPath for processing

    if nargin < 3
        selectedPath = '';
    end

    % make a local copy of options so we can override folderPath if selectedPath supplied
    opts = options;
    if ~isempty(selectedPath)
        % assume the Helper.processRawData uses opts.folderPath (or opts.selectedPath) — override folderPath to selectedPath
        opts.folderPath = selectedPath;
    end

    [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, scan_parameter_names, scan_reference_values] = Helper.processRawData(opts);

    if opts.SAVE_TO_WORKSPACE
        fprintf('\n[INFO] Completed computing OD images. Stored in workspace for reuse.\n');
    else
        fprintf('\n[INFO] Completed computing OD images. Stored on disk for reuse.\n');
    end
end

function fullodimagesFolder = createFullODImagesFolderPath(parentFolder, dataSourcesStruct)
    runID = sprintf('%s_%s_Run%04d', ...
                    dataSourcesStruct{1}.sequence, ...
                    strrep(dataSourcesStruct{1}.date,'/','-'), ...
                    dataSourcesStruct{1}.runs);
    fullodimagesFolder = fullfile(parentFolder, ['FullODImages_' runID]);
end

function [mat_files, raw_scan_parameter_names, raw_scan_parameter_values, raw_file_list, nFiles] = prepareFromOnDiskData(folder)
    mat_files                 = dir(fullfile(folder,'*.mat'));
    mat_files                 = mat_files(~strcmp({mat_files.name},'metadata.mat')); % exclude metadata
    nFiles                    = numel(mat_files);
    raw_scan_parameter_names  = cell(1,nFiles);
    raw_scan_parameter_values = [];
    raw_file_list             = strings(1,nFiles);
end