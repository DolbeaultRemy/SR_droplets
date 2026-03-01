function [full_od_imgs, full_bkg_imgs, raw_scan_parameter_values, raw_file_list, scan_parameter_names, scan_reference_values] = processRawData(options)
%% processRawData
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Reads HDF5 files, computes OD images, supports disk-backed storage in blocks.
%
% Notes:
%   Optional notes, references.

    fprintf('\n[INFO] Processing raw data files at %s ...\n', options.folderPath);

    groupList = ["/images/ODT_1_Axis_Camera/in_situ_absorption", ...
                 "/images/ODT_2_Axis_Camera/in_situ_absorption", ...
                 "/images/Horizontal_Axis_Camera/in_situ_absorption", ...
                 "/images/Vertical_Axis_Camera/in_situ_absorption"];

    % --- Validate camera index ---
    if options.cam < 1 || options.cam > numel(groupList)
        error('Invalid camera index: %d', options.cam);
    end

    files = dir(fullfile(options.folderPath, '*.h5'));
    nFiles = numel(files);
    if nFiles == 0
        error('No HDF5 files found in %s', options.folderPath);
    end

    % Determine image size
    testFile = fullfile(files(1).folder, files(1).name);
    cameraGroup = groupList(options.cam);
    try
        info = h5info(testFile, cameraGroup);      % Check if group exists
    catch
        error('Group "%s" not found in file "%s". Aborting.', cameraGroup, testFile);
    end
    
    datasetNames = {info.Datasets.Name};
    if ~ismember('atoms', datasetNames)
        error('Dataset "%s/atoms" not found in file "%s". Aborting.', cameraGroup, testFile);
    end
    % If we reach here, the dataset exists
    atm_test = double(imrotate(h5read(testFile, append(cameraGroup, "/atoms")), ...
                              options.angle, 'bilinear', 'crop'));
    [ny, nx] = size(atm_test);

    full_od_imgs              = [];
    full_bkg_imgs             = [];
    raw_scan_parameter_values = [];
    raw_file_list             = string(zeros(1,nFiles)); % always string array

    if options.SAVE_TO_WORKSPACE
        fprintf('\n[INFO] Creating in-memory arrays of raw data...\n');
        full_od_imgs  = nan(ny, nx, nFiles, 'single');
        full_bkg_imgs = nan(ny, nx, nFiles, 'single');
    else
        % --- Create uniquely identified full OD image folder ---
        dataSource = makeDataSourceStruct(options.folderPath);
        runID = sprintf('%s_%s_Run%04d', ...
                        dataSource{1}.sequence, ...
                        strrep(dataSource{1}.date,'/','-'), ...
                        dataSource{1}.runs);

        if isfield(options, 'FullODImagesFolder') && ...
           isfolder(options.FullODImagesFolder) && ...
           ~isempty(options.FullODImagesFolder)
            parentFolder = options.FullODImagesFolder;
        else
            parentFolder = options.saveDirectory;
        end

        fullODImageFolder = fullfile(parentFolder, ['FullODImages_' runID]);
        if ~exist(fullODImageFolder,'dir'), mkdir(fullODImageFolder); end
        fprintf('\n[INFO] Creating folder of full OD images on disk: %s\n', fullODImageFolder);

        % --- Save metadata for this run ---
        metadata.options   = options;
        metadata.timestamp = datetime;
        metadata.runID     = runID;
        metadata.imageSize = [ny, nx];
        metadata.fileList  = string(arrayfun(@(f) fullfile(f.folder, f.name), files, 'UniformOutput', false));
        save(fullfile(fullODImageFolder,'metadata.mat'),'metadata','-v7.3');
    end

    % --- Prepare file names ---
    fullFileNames = string(arrayfun(@(f) fullfile(f.folder, f.name), files, 'UniformOutput', false));

    % --- Use specified scan parameter, auto-detect if not specified  ---
    if isfield(options,'scan_parameter') && ~isempty(options.scan_parameter)
        fprintf('\n[INFO] Using user-specified scan parameter(s): ');
        if iscell(options.scan_parameter)
            fprintf('%s\n', strjoin(options.scan_parameter, ', '));
            scan_parameter_names = options.scan_parameter;
            nParams = numel(scan_parameter_names);
        else
            fprintf('%s\n', options.scan_parameter);
            scan_parameter_names = options.scan_parameter;
            nParams = 1;
        end
    else
        [scan_parameter_names_unfiltered, nParams] = detectScanParametersFromFiles(fullFileNames); 
        if isfield(options,'ignore_scan_parameter') && ~isempty(options.ignore_scan_parameter)
            fprintf('\n[INFO] Ignoring the following scan parameter(s): %s\n', strjoin(options.ignore_scan_parameter, ', '));
            scan_parameter_names = scan_parameter_names_unfiltered(~ismember(scan_parameter_names_unfiltered, options.ignore_scan_parameter));
        else
            scan_parameter_names = scan_parameter_names_unfiltered;
        end
    end

    if ~exist('nParams','var')
        nParams = 1;
    end

    % --- Preallocate temp scan values for parfor ---
    temp_scan_vals = cell(nFiles,1);

    % --- Check for Parallel Computing Toolbox ---
    useParallel = license('test','Distrib_Computing_Toolbox') && ~options.SAVE_TO_WORKSPACE;

    if useParallel
        fprintf('\n[INFO] Parallel Computing Toolbox detected. Using parallelization for raw data processing...\n');
        raw_file_list = fullFileNames;

        parfor k = 1:nFiles
            [od_img, bkg_img, val] = readAndComputeOD(fullFileNames(k), options, groupList, ny, nx, scan_parameter_names);
            temp_scan_vals{k} = val(:).';

            if options.SAVE_TO_WORKSPACE
                full_od_imgs(:,:,k)  = single(od_img);
                full_bkg_imgs(:,:,k) = single(bkg_img);
            else
                writeFullODImagesToDisk(fullODImageFolder, od_img, bkg_img, scan_parameter_names, val, fullFileNames(k), k);
            end
        end
    else
        showPB = isfield(options,'showProgressBar') && options.showProgressBar;
        if showPB && options.SAVE_TO_WORKSPACE
            pb = Helper.ProgressBar();
            pb.run('Progress: ');
        end

        for k = 1:nFiles
            [od_img, bkg_img, val] = readAndComputeOD(fullFileNames(k), options, groupList, ny, nx, scan_parameter_names);
            temp_scan_vals{k} = val(:).';

            if options.SAVE_TO_WORKSPACE
                full_od_imgs(:,:,k)  = single(od_img);
                full_bkg_imgs(:,:,k) = single(bkg_img);
            else
                writeFullODImagesToDisk(fullODImageFolder, od_img, bkg_img, scan_parameter_names, val, fullFileNames(k), k);
            end
            raw_file_list(k) = fullFileNames(k);

            if showPB && options.SAVE_TO_WORKSPACE
                progressPercent = round(k/nFiles*100);
                pb.run(progressPercent);
            end
        end

        if showPB && options.SAVE_TO_WORKSPACE
            pb.run(' Done!');
        end
    end
    
    % --- Convert cell array to numeric matrix after parfor ---
    if nParams == 1
        raw_scan_parameter_values = cellfun(@double, temp_scan_vals);                          % enforce double here
        raw_scan_parameter_values = reshape(raw_scan_parameter_values, 1, []);                 % row vector
    else
        raw_scan_parameter_values = cellfun(@double, temp_scan_vals, 'UniformOutput', false);
        raw_scan_parameter_values = vertcat(raw_scan_parameter_values{:});                     % rows = files, cols = parameters
    end
    
    % --- Determine scan reference values ---
    if (~isfield(options,'scan_reference_values') || isempty(options.scan_reference_values))
        if isvector(raw_scan_parameter_values)
            % Single parameter case
            scan_reference_values = unique(raw_scan_parameter_values(:), 'stable');
            scan_reference_values = scan_reference_values(:).'; % row vector
        else
            % Multi-parameter case: unique rows
            scan_reference_values = unique(raw_scan_parameter_values, 'rows', 'stable');
        end
    else
        % Ensure anything coming from options is also double **once, here**
        scan_reference_values = double(options.scan_reference_values);
    end
end

%% --- Local helper functions ---
function [od_img, bkg_img, val] = readAndComputeOD(fullFileName, options, groupList, ny, nx, scanParamNames)
    try
        atm_img  = double(imrotate(h5read(fullFileName, append(groupList(options.cam), "/atoms")), options.angle, 'bilinear', 'crop'));
        bkg_img  = double(imrotate(h5read(fullFileName, append(groupList(options.cam), "/background")), options.angle, 'bilinear', 'crop'));
        dark_img = double(imrotate(h5read(fullFileName, append(groupList(options.cam), "/dark")), options.angle, 'bilinear', 'crop'));
        od_img   = Helper.calculateODImage(atm_img, bkg_img, dark_img, options.ImagingMode, options.PulseDuration);
    catch
        warning('\nMissing data in %s, storing NaNs.', fullFileName);
        od_img  = nan(ny, nx);
        bkg_img = nan(ny, nx);
    end

    % --- Read scan parameter(s) for this file ---
    try
        if iscell(scanParamNames)
            val = NaN(1,numel(scanParamNames));
            for j = 1:numel(scanParamNames)
                val(j) = h5readatt(fullFileName, '/globals', scanParamNames{j});
            end
        else
            val = h5readatt(fullFileName, '/globals', scanParamNames);
        end
    catch
        val = NaN;
    end
end

function [scanParamNames, nParams] = detectScanParametersFromFiles(fileNames, minFilesToCheck)
% Detect scan parameter(s) by checking which numeric attributes vary across multiple files
%
% Inputs:
%   fileNames        - string array or cell array of full HDF5 file paths
%   minFilesToCheck  - minimum number of files to compare (default: 3)
%
% Outputs:
%   scanParamNames   - char array (single) or cell array of char arrays (multiple)
%   nParams          - number of detected scan parameters

    if nargin < 2
        minFilesToCheck = 3;
    end

    % Ensure fileNames is a string array
    if iscell(fileNames)
        fileNames = string(fileNames);
    end

    nFiles = numel(fileNames);
    nCheck = min(minFilesToCheck, nFiles);
    fprintf('\n[INFO] Detecting scan parameter(s)...\n');

    % Read attribute names from first file
    info = h5info(fileNames(1));
    if any(strcmp({info.Groups.Name}, '/globals'))
        globalsGroup = info.Groups(strcmp({info.Groups.Name}, '/globals'));
    else
        warning('\n[WARNING] /globals group not found in first file.');
        scanParamNames = NaN;
        nParams = 0;
        return;
    end

    attrNames = string({globalsGroup.Attributes.Name});
    numericAttrNames = strings(0,1);

    % Identify numeric attributes
    for j = 1:numel(attrNames)
        if attrNames(j) == "runs"
            continue;
        end
        val = h5readatt(fileNames(1), '/globals', attrNames(j));
        if isnumeric(val)
            numericAttrNames(end+1) = attrNames(j); %#ok<AGROW>
        end
    end

    if isempty(numericAttrNames)
        fprintf('\n[WARNING] No numeric attributes in /globals.\n');
        scanParamNames = NaN;
        nParams = 0;
        return;
    end

    % --- Check which numeric attributes vary across first few files ---
    varyingParams = strings(0,1);
    for j = 1:numel(numericAttrNames)
        attrName = numericAttrNames(j);
        values = NaN(1, nCheck);
        for k = 1:nCheck
            try
                val = h5readatt(fileNames(k), '/globals', attrName);
                if isnumeric(val)
                    values(k) = val(1); % use first element if array
                end
            catch
                values(k) = NaN;
            end
        end
        % Ignore NaNs when checking for variation
        if numel(unique(values(~isnan(values)))) > 1
            varyingParams(end+1) = attrName; %#ok<AGROW>
        end
    end

    % --- Return result ---
    nParams = numel(varyingParams);
    if nParams == 0
        fprintf('\n[INFO] No varying scan parameters detected.\n');
        scanParamNames = NaN;
    elseif nParams == 1
        scanParamNames = char(varyingParams); % single char array
        fprintf('\n[INFO] Single scan parameter detected: %s\n', scanParamNames);
    else
        scanParamNames = cellstr(varyingParams); % cell array of char arrays
        fprintf('\n[INFO] Multiple scan parameters detected: %s\n', strjoin(varyingParams, ', '));
    end
end

function writeFullODImagesToDisk(fullODImageFolder, od_img, bkg_img, param, val, file_name, idx)
% Writes OD/BKG image + scan parameter(s) to a MAT file
    matFilePath = fullfile(fullODImageFolder, sprintf('Image_%04d.mat', idx));
    OD          = single(od_img);
    BKG         = single(bkg_img);
    File        = string(file_name);
    Scan_Param  = param;
    Scan_Val    = single(val);
    save(matFilePath, 'OD','BKG', 'Scan_Param', 'Scan_Val','File','-v7.3');
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
