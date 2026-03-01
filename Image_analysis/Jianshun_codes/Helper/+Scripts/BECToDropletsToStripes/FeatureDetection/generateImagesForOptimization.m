%% generateImagesForOptimization.m
% Converts a cell array of images into PNG files for Bayesian optimization
% Launches Image Segmenter for manually labeling ground truth masks

% ----------------- USER INPUTS -----------------
outputFolder = 'OptimizationImages'; % folder to save PNGs
maskFolder   = 'OptimizationMasks';  % folder to save masks
startIdx     = 1;           % first image to export
endIdx       = 10;           % last image to export (set = 1 for a single test image)
createMasks  = true;        % set true to generate masks manually
% -----------------------------------------------

% Create folders if they don't exist
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
if createMasks && ~exist(maskFolder,'dir')
    mkdir(maskFolder);
end

% Loop over selected images
for idx = startIdx:endIdx
    img = od_imgs{idx};
    
    % Convert to double if not already
    if ~isa(img,'double')
        img = im2double(img);
    end
    
    % Normalize to [0,1] for PNG
    img = img - min(img(:));
    img = img / max(img(:));
    
    % Save image
    imgName = sprintf('image_%03d.png', idx);
    imwrite(img, fullfile(outputFolder, imgName));
    
   % Optional: create ground truth mask using Image Segmenter
    if createMasks
        fprintf('\n[INFO] Segment image %d/%d using Image Segmenter...\n', idx, endIdx);
    
        % Launch Image Segmenter with the image
        imageSegmenter(img);
    
        % --- Instructions for user ---
        % 1. Segment the regions you want as ground truth.
        % 2. Click "Export" → "Export to Workspace", export the binary mask
        % 3. Close the Image Segmenter when done.
        % ------------------------------------------------
    
        % Pause script until user confirms they exported the mask
        input('\n[INFO] After exporting mask from Image Segmenter, press Enter to continue...','s');
        
        % Look for any BW variables in base workspace
        bwVars = evalin('base', 'who(''BW*'')');
        if ~isempty(bwVars)
            % Pick the last one (e.g. BW, BW1, BW2, ...)
            latestBW = bwVars{end};
            gtMask = evalin('base', latestBW);
            
            % Save mask
            maskName = sprintf('mask_%03d.png', idx);
            imwrite(gtMask, fullfile(maskFolder, maskName));
            fprintf('\n[INFO] Mask saved: %s (from %s)\n', maskName, latestBW);
            
            % Clear masks from workspace to avoid confusion
            for v = 1:numel(bwVars)
                evalin('base', sprintf('clear %s', bwVars{v}));
            end
        else
            warning('No mask found. Skipping mask save for image %d.', idx);
        end
    end
end

fprintf('\n[INFO] Image and mask export complete.\n');

