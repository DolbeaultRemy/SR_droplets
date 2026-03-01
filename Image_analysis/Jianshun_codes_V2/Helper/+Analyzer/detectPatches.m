function [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart] = detectPatches(img, params)
%% detectPatches
% Author:       Karthik
% Date:         2025-09-12
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
    
    %% --- Step 1: Background subtraction & cloud mask ---
    % Morphological opening estimates slowly-varying background
    seRadius = max(3, round(min(size(img)) * params.backgroundDiskFraction));
    backgroundEstimate = imopen(img, strel('disk', seRadius));
    
    % Subtract background and clamp negatives to zero
    imgCorrected = img - backgroundEstimate;
    imgCorrected(imgCorrected < 0) = 0;
    
    %% --- Step 2: Cloud segmentation ---
    % Smooth image to remove high-frequency noise
    imgSmoothed = imgaussfilt(imgCorrected, round(min(size(img))/25));
    
    % Threshold using Otsu to create binary cloud mask
    cloudMask = imbinarize(imgSmoothed, graythresh(imgSmoothed));
    
    % Close small gaps and fill holes
    cloudMask = imclose(cloudMask, strel('disk', round(seRadius/4)));
    cloudMask = imfill(cloudMask,'holes');
    
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
    
    patchProps = detectPatchesCore(imgDenoised, opts);
    
    if ~isempty(patchProps)
        patchCentroidsLocal  = cat(1, patchProps.Centroid);
        patchCentroidsGlobal = patchCentroidsLocal + [xStart-1, yStart-1];
    else
        patchCentroidsGlobal = [];
    end
end

%% --- Helper function ---
function patchProps = detectPatchesCore(I, opts)
% Detect lattice patches (blobs/stripes) using Difference-of-Gaussians,
% thresholding, and connected component analysis.
% Returns a struct array with fields:
%   Centroid, Area, Orientation, MajorAxisLength, MinorAxisLength
    
    % Step 1: Difference-of-Gaussians
    G1          = imgaussfilt(I, opts.sigmaSmall);
    G2          = imgaussfilt(I, opts.sigmaLarge);
    dogResponse = mat2gray(G1 - G2);
    
    % Step 2: Adaptive threshold
    T           = adaptthresh(dogResponse, opts.adaptiveSensitivity, 'NeighborhoodSize', opts.adaptiveNeighborhoodSize, 'Statistic','gaussian');
    binaryMask  = imbinarize(dogResponse, T);
    
    % Step 3: Remove small specks, close small gaps
    binaryMask  = bwareaopen(binaryMask, 5);
    binaryMask  = imclose(binaryMask, strel('disk',1));
    
    % Step 4: Connected components
    CC = bwconncomp(binaryMask);
    if CC.NumObjects == 0
        patchProps = struct('Centroid',[],'Area',[],'Orientation',[], ...
                            'MajorAxisLength',[],'MinorAxisLength',[]);
        return;
    end
    
    % Step 5: Extract patch properties
    patchPropsAll = regionprops(CC, dogResponse, 'Centroid','Area','Orientation','MajorAxisLength','MinorAxisLength');
    
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
end
