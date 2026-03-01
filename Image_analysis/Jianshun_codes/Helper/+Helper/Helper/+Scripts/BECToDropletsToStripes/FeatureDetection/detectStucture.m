function [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart, binaryMask, CC] = detectStucture(img, params, folderPath, fileNum)
%% detectPatches
% Author:       Jianshun
% Date:         2025-09-27
% Version:      1.0
%
% Description:
%   Detect lattice patches (blobs/stripes) in a single OD image.
%   Performs background subtraction, cloud segmentation, cropping,
%   denoising, Difference-of-Gaussians filtering, patch detection, and plotting.
%
% Inputs:
%   img    - 2D OD image (double or converted internally)
%   params - struct of user-tunable parameters:
%       backgroundDiskFraction - fraction of image size for morphological opening
%       boundingBoxPadding     - pixels of padding around cloud bounding box
%       dogGaussianSmallSigma  - sigma for small Gaussian in DoG
%       dogGaussianLargeSigma  - sigma for large Gaussian in DoG
%       minPeakProminence      - min DoG response for thresholding
%       minPeakFraction        - fraction of max DoG response for adaptive threshold
%       subpixelWindowRadius   - radius for potential subpixel refinement (not used here)
%       minimumPatchArea       - minimum patch area to keep
%       pixelSize              - meters/pixel
%       magnification          - imaging system magnification
%   hAx    - axes handle for plotting
%
% Outputs:
%   patchProps           - struct array of detected patches
%   patchCentroidsGlobal - Nx2 array of patch centroids in image coordinates
%   imgCropped           - cropped, background-subtracted image used for detection
%   xAxis, yAxis         - physical axes in microns
%
% Notes:
%   Optional notes, references.

    if ~isa(img,'double')
        img = im2double(img);
    end
    [Ny, Nx] = size(img);

    binaryMask = [];

    figure(1)
    imagesc(img, [0, 4])
    colormap jet
    
    save_path = folderPath + "AnalysisPlot\";
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    % filename = 'RawImg' + str(fileNum) + '.jpg';
    % full_path = fullfile(save_path, filename);
    % % 导出图像
    % print(gcf, full_path, '-djpeg', '-r300'); 
    
    %% --- Step 1: Background subtraction & cloud mask ---
    % Morphological opening estimates slowly-varying background
    seRadius = max(3, round(min(size(img)) * params.backgroundDiskFraction));
    backgroundEstimate = imopen(img, strel('disk', seRadius));
    
    % Subtract background and clamp negatives to zero
    imgCorrected = img - backgroundEstimate;
    imgCorrected(imgCorrected < 0) = 0;

    figure(2)
    imagesc(imgCorrected, [0, 4])
    colormap jet
    
    %% --- Step 2: Cloud segmentation ---
    % Smooth image to remove high-frequency noise
    imgSmoothed = imgaussfilt(imgCorrected, round(min(size(img))/100));

    figure(3)
    imagesc(imgSmoothed, [0, 4])
    colormap jet
    
    % Threshold using Otsu to create binary cloud mask
    cloudMask = imbinarize(imgSmoothed, graythresh(imgSmoothed));

    figure(4)
    imagesc(cloudMask, [0, 4])
    colormap jet
    
    % Close small gaps and fill holes
    % cloudMask = imclose(cloudMask, strel('disk', round(seRadius/0.5)));
    % cloudMask = imfill(cloudMask,'holes');
    
    %% --- Step 3: Largest connected region & crop ---
    CC = bwconncomp(cloudMask);
    if CC.NumObjects > 0
        stats = regionprops(CC,'Area','BoundingBox');
        [~, idxMax] = max([stats.Area]);
        bb = round(stats(idxMax).BoundingBox);
    
        % Crop with padding, stay within image
        xStart = max(1, bb(1)-params.boundingBoxPadding);
        yStart = max(1, bb(2)-params.boundingBoxPadding);
        xEnd   = min(Nx, bb(1)+bb(3)-1 + params.boundingBoxPadding);
        yEnd   = min(Ny, bb(2)+bb(4)-1 + params.boundingBoxPadding);
    else
        % Fallback: full image
        xStart=1; yStart=1; xEnd=Nx; yEnd=Ny;
    end
    imgCropped = imgCorrected(yStart:yEnd, xStart:xEnd);

    figure(5)
    imagesc(imgCropped, [0, 4])
    colormap jet
    % filename = "RawImg" + fileNum + ".jpg";
    % full_path = fullfile(save_path, filename);
    % % 导出图像
    % print(gcf, full_path, '-djpeg', '-r300'); 
    
    %% --- Step 4: Denoising ---
    % Light Gaussian blur to reduce high-frequency noise
    imgDenoised = imgaussfilt(imgCropped, 0.8);
    
    %% --- Step 5: Detect lattice patches ---
    opts = struct('sigmaSmall', params.dogGaussianSmallSigma, ...
                  'sigmaLarge', params.dogGaussianLargeSigma, ...
                  'adaptiveSensitivity', params.adaptiveSensitivity, ...
                  'adaptiveNeighborhoodSize', params.adaptiveNeighborhoodSize, ...
                  'minPatchArea', params.minimumPatchArea, ...
                  'minPeakFraction', params.minPeakFraction);
    
    [patchProps, binaryMask, CC] = detectPatchesCore(imgDenoised, opts, imgCropped, save_path, fileNum);
    
    if ~isempty(patchProps)
        patchCentroidsLocal  = cat(1, patchProps.Centroid);
        patchCentroidsGlobal = patchCentroidsLocal + [xStart-1, yStart-1];
    else
        patchCentroidsGlobal = [];
    end
end

%% --- Helper function ---
function [patchProps, binaryMask, CC] = detectPatchesCore(I, opts, img, folderPath, fileNum)
% Detect lattice patches (blobs/stripes) using Difference-of-Gaussians,
% thresholding, and connected component analysis.
% Returns a struct array with fields:
%   Centroid, Area, Orientation, MajorAxisLength, MinorAxisLength

    % Step 0: Smooth the imaging
    imgSmoothed = imgaussfilt(I, 0.1);
    figure(8)
    imagesc(imgSmoothed, [0, 4])
    colormap jet

    temp = I;
    I = imgSmoothed;

    
    % Step 1: Difference-of-Gaussians
    G1          = imgaussfilt(I, opts.sigmaSmall); %0.5
    G2          = imgaussfilt(I, opts.sigmaLarge); %4
    dogResponse = mat2gray(G1 - G2);

    figure(6)
    imagesc(dogResponse, [0, 1])
    colormap jet
    
    % Step 2: Adaptive threshold
    T           = adaptthresh(dogResponse, opts.adaptiveSensitivity, 'NeighborhoodSize', opts.adaptiveNeighborhoodSize, 'Statistic','gaussian');
    binaryMask  = imbinarize(dogResponse, T);

    figure(7)
    imagesc(binaryMask, [0, 1])
    colormap jet 
    
    % Step 3: Remove small specks, close small gaps
    binaryMask  = bwareaopen(binaryMask, 5);
    binaryMask  = imclose(binaryMask, strel('disk',4));

    % figure(9)
    % imagesc(binaryMask, [0, 1])
    % colormap jet
    % filename = "Mask" + fileNum + ".jpg";
    % full_path = fullfile(folderPath, filename);
    % % 导出图像
    % print(gcf, full_path, '-djpeg', '-r300');
    
    % Step 4: Connected components
    CC = bwconncomp(binaryMask);
    if CC.NumObjects == 0
        patchProps = struct('Centroid',[],'Area',[],'Orientation',[], ...
                            'MajorAxisLength',[],'MinorAxisLength',[], ...
                            'PixelIdxList', []);
        return;
    end

    figure(10)
    imagesc(img, [0, 3])
    colormap jet
    hold on;
    
    % Step 5: Extract patch properties
    patchPropsAll = regionprops(CC, dogResponse, 'Centroid','Area','Orientation','MajorAxisLength','MinorAxisLength','PixelIdxList');
    patchProps.Connectivity = CC.Connectivity;
    patchProps.ImageSize = CC.ImageSize; 
    patchProps.NumObjects = CC.NumObjects;
    % patchPropsAll.PixelIdxList = CC.PixelIdxList;
    
    % Step 6: Filter by minimum area and adaptive per-patch peak fraction
    maxDog = max(dogResponse(:));                    % Global maximum DoG response
    keepIdx = false(1,numel(patchPropsAll));
    for k = 1:numel(patchPropsAll)
        patchPixels = CC.PixelIdxList{k};
        patchMax = max(dogResponse(patchPixels));    % Peak DoG in this patch
        if patchPropsAll(k).Area >= opts.minPatchArea && patchMax >= opts.minPeakFraction*maxDog
            keepIdx(k) = true;
        end
    end
    patchProps = patchPropsAll(keepIdx);

    bwImage = false(CC.ImageSize);
    for i = 1:CC.NumObjects
        if keepIdx(i)
            bwImage(CC.PixelIdxList{i}) = true;
        end
    end

    CC = bwconncomp(bwImage);
    binaryMask = bwImage;
    
    boundaries = bwboundaries(bwImage);
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1.5);
    end
    
    title(sprintf('Find %d Areas', CC.NumObjects));
    legend('boundary');
    hold off;

    % save_path = folderPath + "Detect\Update\";
    % if ~exist(save_path, 'dir')
    %     mkdir(save_path);
    % end
    % filename = "Detect" + fileNum + ".jpg";
    % full_path = fullfile(save_path, filename);
    % % 导出图像
    % print(gcf, full_path, '-djpeg', '-r300');

end
