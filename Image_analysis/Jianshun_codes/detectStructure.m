function [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart, binaryMask, CC] = ...
    detectStructure(img, params, folderPath, fileName)
% DETECTSTRUCTURE Detect lattice patches (blobs/stripes) in an OD image.
% Inputs:
%   img        - 2D optical density image (double)
%   params     - struct with detection parameters
%   folderPath - string, path to save intermediate plots
%   fileName   - string, base name for saved figures
% Outputs:
%   patchProps           - struct array of detected patches
%   patchCentroidsGlobal - Nx2 centroids in original image coordinates
%   imgCropped           - cropped background-subtracted image
%   xStart, yStart       - cropping offsets
%   binaryMask           - binary mask of detected patches
%   CC                   - connected components of binaryMask

    if ~isa(img, 'double')
        img = im2double(img);
    end
    [Ny, Nx] = size(img);

    % Display raw image
    figure(1); clf;
    imagesc(img, [0, 4]);
    colormap jet;
    title('Raw OD Image');
    colorbar;

    %% Step 1: Background subtraction using morphological opening
    seRadius = max(3, round(min(size(img)) * params.backgroundDiskFraction));
    backgroundEstimate = imopen(img, strel('disk', seRadius));
    imgCorrected = img - backgroundEstimate;
    imgCorrected(imgCorrected < 0) = 0;

    figure(2); clf;
    imagesc(imgCorrected, [0, 4]);
    colormap jet;
    title('Background Subtracted');

    %% Step 2: Cloud segmentation (smooth and threshold)
    smoothingSigma = round(min(size(img)) / 100);
    imgSmoothed = imgaussfilt(imgCorrected, smoothingSigma);
    figure(3); clf;
    imagesc(imgSmoothed, [0, 4]);
    colormap jet;
    title('Smoothed Image');

    cloudMask = imbinarize(imgSmoothed, graythresh(imgSmoothed));
    figure(4); clf;
    imagesc(cloudMask);
    colormap gray;
    title('Initial Cloud Mask');

    %% Step 3: Crop to largest connected component
    CC_cloud = bwconncomp(cloudMask);
    if CC_cloud.NumObjects > 0
        stats = regionprops(CC_cloud, 'Area', 'BoundingBox');
        [~, idxMax] = max([stats.Area]);
        bb = round(stats(idxMax).BoundingBox);
        xStart = max(1, bb(1) - params.boundingBoxPadding);
        yStart = max(1, bb(2) - params.boundingBoxPadding);
        xEnd   = min(Nx, bb(1) + bb(3) - 1 + params.boundingBoxPadding);
        yEnd   = min(Ny, bb(2) + bb(4) - 1 + params.boundingBoxPadding);
    else
        xStart = 1; yStart = 1; xEnd = Nx; yEnd = Ny;
    end
    imgCropped = imgCorrected(yStart:yEnd, xStart:xEnd);

    figure(5); clf;
    imagesc(imgCropped, [0, 4]);
    colormap jet;
    title('Cropped Cloud Region');

    %% Step 4: Denoising
    imgDenoised = imgaussfilt(imgCropped, 0.8);

    %% Step 5: Detect patches using Difference-of-Gaussians
    opts = struct('sigmaSmall', params.dogGaussianSmallSigma, ...
                  'sigmaLarge', params.dogGaussianLargeSigma, ...
                  'adaptiveSensitivity', params.adaptiveSensitivity, ...
                  'adaptiveNeighborhoodSize', params.adaptiveNeighborhoodSize, ...
                  'minPatchArea', params.minimumPatchArea, ...
                  'minPeakFraction', params.minPeakFraction);

    [patchProps, binaryMask, CC] = detectPatchesCore(imgDenoised, opts, imgCropped, folderPath, fileName);

    if ~isempty(patchProps)
        patchCentroidsLocal = cat(1, patchProps.Centroid);
        patchCentroidsGlobal = patchCentroidsLocal + [xStart-1, yStart-1];
    else
        patchCentroidsGlobal = [];
    end
end

%% ------------------------------------------------------------------------
function [patchProps, binaryMask, CC] = detectPatchesCore(I, opts, imgOriginal, folderPath, fileName)
% detectPatchesCore: Core patch detection using DoG and adaptive threshold.
    % Slight smoothing
    imgSmoothed = imgaussfilt(I, 0.1);
    figure(8); clf;
    imagesc(imgSmoothed, [0, 4]);
    colormap jet;
    title('Smoothed for DoG');

    I = imgSmoothed; % use smoothed for DoG

    % Difference-of-Gaussians
    G1 = imgaussfilt(I, opts.sigmaSmall);
    G2 = imgaussfilt(I, opts.sigmaLarge);
    dogResponse = mat2gray(G1 - G2);

    figure(6); clf;
    imagesc(dogResponse, [0, 1]);
    colormap jet;
    title('DoG Response');

    % Adaptive threshold
    T = adaptthresh(dogResponse, opts.adaptiveSensitivity, ...
        'NeighborhoodSize', opts.adaptiveNeighborhoodSize, 'Statistic', 'gaussian');
    binaryMask = imbinarize(dogResponse, T);

    figure(7); clf;
    imagesc(binaryMask);
    colormap gray;
    title('Initial Binary Mask');

    % Clean mask: remove small specs, close gaps
    binaryMask = bwareaopen(binaryMask, 5);
    binaryMask = imclose(binaryMask, strel('disk', 4));

    % Connected components
    CC = bwconncomp(binaryMask);
    if CC.NumObjects == 0
        patchProps = struct('Centroid', [], 'Area', [], 'Orientation', [], ...
                            'MajorAxisLength', [], 'MinorAxisLength', [], ...
                            'PixelIdxList', []);
        return;
    end

    % Extract properties
    patchPropsAll = regionprops(CC, dogResponse, 'Centroid', 'Area', ...
        'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'PixelIdxList');

    % Filter by area and peak response
    maxDog = max(dogResponse(:));
    keepIdx = false(1, numel(patchPropsAll));
    for k = 1:numel(patchPropsAll)
        patchPixels = CC.PixelIdxList{k};
        patchMax = max(dogResponse(patchPixels));
        if patchPropsAll(k).Area >= opts.minPatchArea && patchMax >= opts.minPeakFraction * maxDog
            keepIdx(k) = true;
        end
    end
    patchProps = patchPropsAll(keepIdx);

    % Create new binary mask and CC from kept patches
    bwImage = false(CC.ImageSize);
    for i = 1:CC.NumObjects
        if keepIdx(i)
            bwImage(CC.PixelIdxList{i}) = true;
        end
    end
    CC = bwconncomp(bwImage);
    binaryMask = bwImage;

    % Plot boundaries on original cropped image
    figure(10); clf;
    imagesc(imgOriginal, [0, 3]);
    colormap jet;
    hold on;
    boundaries = bwboundaries(binaryMask);
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1.5);
    end
    title(sprintf('Detected %d Patches', CC.NumObjects));
    legend('boundary');
    hold off;

    % Optionally save figure
    % saveFigure(folderPath, fileName, 'Detect');
end

% function saveFigure(folderPath, fileName, subfolder)
%     savePath = fullfile(folderPath, 'AnalysisPlot', subfolder);
%     if ~exist(savePath, 'dir'); mkdir(savePath); end
%     print(gcf, fullfile(savePath, ['Detect_' fileName '.jpg']), '-djpeg', '-r300');
% end