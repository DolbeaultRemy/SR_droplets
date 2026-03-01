function [SAVE_TO_WORKSPACE, runMemoryGB] = estimateDatasetMemory(dataSources, options)
%% estimateDatasetMemory
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Estimate per-run memory and decide whether to save dataset to workspace
%   Supports both raw data folders and a preselected run folder.
%
% Notes:
%   Optional notes, references.

    % --- Measured memory per image (bytes) ---
    bytesPerFullImage    = 37.75 * 1e6;  % full OD image
    bytesPerCroppedImage = 0.16  * 1e6;  % cropped OD image

    % --- Check available RAM on Windows ---
    if ispc
        [~, sys] = memory;
        availableRAM = sys.PhysicalMemory.Available;
    else
        availableRAM = 8e9; % fallback: 8 GB if not Windows
    end

    SAVE_TO_WORKSPACE = true; % default, may change per run
    runMemoryGB = [];         % store per-run memory

    % --- Case 1: selected path exists ---
    if isfield(options, 'selectedPath') && isfolder(options.selectedPath)
        runFolder = options.selectedPath;
        files = dir(fullfile(runFolder, '*.h5'));
        nFiles = numel(files);

        if nFiles > 0
            runBytes = nFiles * (bytesPerFullImage + bytesPerCroppedImage);
            runMemoryGB(end+1,1) = runBytes/1e9;

            if runBytes > 0.75 * availableRAM
                SAVE_TO_WORKSPACE = false;
                fprintf('\n[INFO] Selected run %s estimated size %.2f GB exceeds 75%% of available RAM. Will save images to disk if not already available.\n', ...
                        runFolder, runBytes/1e9);
            else
                fprintf('\n[INFO] Selected run %s estimated size %.2f GB fits in memory.\n', ...
                        runFolder, runBytes/1e9);
            end
        end

    % --- Case 2: fallback to raw data folders ---
    else
        for i_ds = 1:numel(dataSources)
            ds = dataSources{i_ds};
            if isfield(ds, 'baseFolder') && ~isempty(ds.baseFolder)
                baseFolder = fullfile(ds.baseFolder, ds.sequence, ds.date);
            else
                baseFolder = fullfile(options.baseDataFolder, ds.sequence, ds.date);
            end

            for j_run = 1:numel(ds.runs)
                runItem = ds.runs(j_run);
                runID = sprintf('%04d', runItem);
                runFolder = fullfile(baseFolder, runID);

                if isfolder(runFolder)
                    files = dir(fullfile(runFolder, '*.h5'));
                    nFiles = numel(files);
                    if nFiles == 0
                        continue;
                    end

                    runBytes = nFiles * (bytesPerFullImage + bytesPerCroppedImage);
                    runMemoryGB(end+1,1) = runBytes/1e9;

                    if runBytes > 0.75 * availableRAM
                        SAVE_TO_WORKSPACE = false;
                        fprintf('\n[INFO] Run %s/%s estimated size %.2f GB exceeds 75%% of available RAM. Will save images to disk if not already available.\n', ...
                            ds.sequence, runID, runBytes/1e9);
                    else
                        fprintf('\n[INFO] Run %s/%s estimated size %.2f GB fits in memory.\n', ...
                            ds.sequence, runID, runBytes/1e9);
                    end
                end
            end
        end
    end
end
