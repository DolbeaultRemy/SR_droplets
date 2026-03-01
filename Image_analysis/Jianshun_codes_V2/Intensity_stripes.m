addpath(genpath("."));

folderPath = "Y:\StructuralPhaseTransition\2026\02\11\";
ROIcenter = [1500 1950];
measurementName = "Stripes_90°";
run = "0001";
cam = 5;
groupList    = ["/images/MOT_3D_Camera/in_situ_absorption", "/images/ODT_1_Axis_Camera/in_situ_absorption", "/images/ODT_2_Axis_Camera/in_situ_absorption", "/images/Horizontal_Axis_Camera/in_situ_absorption", "/images/Vertical_Axis_Camera/in_situ_absorption"];
span         = [400, 400];
center       = [1500, 1950];
fraction     = [0.1, 0.1];
shotwindow   = 5; %380;
removeFringes = false;


[od_imgs,Date] = run_data_remy(folderPath, run, cam, groupList, span, center, fraction ,shotwindow, removeFringes);

% 
% detectionParams = setDetectionParameters();
% 
% folderPath     = strcat(folderPath, run);    
% filePattern = fullfile(folderPath, '*.h5');
% files       = dir(filePattern);
% 
% 
% %%% Detect structures and create masked images %%%
% % Initialize cell arrays to store results
% masks_detected_structures = {};
% xStart_mask = {};
% yStart_mask = {};
% od_imgs_crop = {};
% masked_images_atoms = {};

% % Process all images
% for idx = 1:length(od_imgs)
%     % Detect structures (your existing code)
%     [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart, binaryMask, CC] = ...
%         detectStructure(od_imgs{1, idx}, detectionParams, folderPath, files(idx).name);
% 
%     % Store detection results
%     masks_detected_structures{end+1} = binaryMask;
%     xStart_mask{end+1} = xStart;
%     yStart_mask{end+1} = yStart;
%     % 
% end

% % Save masked_images to HDF5 file with custom name %%%
% 
% % Set your desired filename here
% customFilename = 'my_masked_images.h5';  % CHANGE THIS TO YOUR DESIRED NAME
% 
% % Create full path
% outputFilename = fullfile(customFilename);
% 
% % Convert to character vector if needed
% outputFilename = char(outputFilename);
% 
% % Delete existing file if it exists
% if exist(outputFilename, 'file') == 2
%     delete(outputFilename);
%     fprintf('Deleted existing file: %s\n', outputFilename);
% end
% 
% % Save each masked image as a separate dataset
% for idx = 1:length(masked_images)
%     % Create dataset name
%     datasetName = sprintf('/image_%04d', idx);
%     datasetName = char(datasetName);
% 
%     % Save the image
%     h5create(outputFilename, datasetName, size(masked_images{idx}), 'Datatype', 'double');
%     h5write(outputFilename, datasetName, masked_images{idx});
% 
%     % Optional: save minimal metadata
%     h5writeatt(outputFilename, datasetName, 'original_index', idx);
% end
% 
% % Save basic info about the dataset
% h5writeatt(outputFilename, '/', 'num_images', length(masked_images));
% h5writeatt(outputFilename, '/', 'description', 'Masked images with non-detected pixels set to 0');
% 
% fprintf('Saved %d masked images to: %s\n', length(masked_images), outputFilename);




