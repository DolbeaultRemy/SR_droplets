function results = extractAndClassifyShapes(imgCropped, patchProps, xStart, yStart, params)
%% extractAndClassifyShapes
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Extracts internal shapes inside DoG-detected patches using intensity
%   thresholding, optional multi-scale DoG fusion, and Canny edge detection.
%   This function is robust to irregular shapes, filaments, and Y-shaped structures.
%
% Inputs:
%   imgCropped    - cropped image (DoG input)
%   patchProps    - struct array from detectPatchesDoG containing patch info
%   xStart, yStart- offsets of the cropped region relative to full image
%   params        - struct with the following fields:
%       shapeMinArea         - minimum area of shapes in pixels
%       shapeCloseRadius     - radius for morphological closing
%       shapeFillHoles       - boolean; fill internal holes in shapes
%       shapeSkeletonize     - boolean; compute skeletons of shapes
%       intensityThreshFraction - fraction of max intensity to threshold
%       edgeSigma            - Gaussian smoothing sigma for Canny
%       edgeThresholdLow     - low threshold fraction for Canny
%       edgeThresholdHigh    - high threshold fraction for Canny
%
% Outputs:
%   results - struct array with fields:
%       BW          - binary mask of detected shapes (local to cropped ROI)
%       labels      - labeled mask of shapes
%       props       - regionprops with classification
%       skeletons   - skeletons of shapes (if enabled)
%       boundaries  - global boundaries in full image coordinates
%
% Notes:
%   Optional notes, references.

    %% --- Return empty if no patches detected ---
    if isempty(patchProps)
        results = struct('BW',{},'labels',{},'props',{},'skeletons',{},'boundaries',{});
        return;
    end

    %% --- Preallocate result structure ---
    results = repmat(struct('BW',[],'labels',[],'props',[],'skeletons',[],'boundaries',[]), numel(patchProps), 1);

    %% --- Loop over each detected patch ---
    for k = 1:numel(patchProps)

        % --- Create local patch mask ---
        % Construct an elliptical mask corresponding to the detected patch
        mask = false(size(imgCropped));
        cx   = patchProps(k).Centroid(1);       % patch center x
        cy   = patchProps(k).Centroid(2);       % patch center y
        a    = patchProps(k).MajorAxisLength/2; % semi-major axis
        b    = patchProps(k).MinorAxisLength/2; % semi-minor axis
        phi  = deg2rad(-patchProps(k).Orientation); % rotation angle in radians

        [xx,yy] = meshgrid(1:size(imgCropped,2), 1:size(imgCropped,1));
        Xr = (xx-cx)*cos(phi) + (yy-cy)*sin(phi);
        Yr = -(xx-cx)*sin(phi) + (yy-cy)*cos(phi);
        mask((Xr/a).^2 + (Yr/b).^2 <= 1) = true;  % points inside ellipse

        % --- Apply mask to ROI ---
        % Only consider pixels inside patch ellipse
        roi = imgCropped;
        roi(~mask) = 0;

        % --- Multi-scale Canny edge detection ---
        % Detect edges to help capture fine structure, filaments, and Y-shaped arms
        edges = edge(roi, 'Canny', [params.edgeThresholdLow, params.edgeThresholdHigh], params.edgeSigma);

        % --- Intensity thresholding ---
        % Retain pixels above a fraction of the maximum intensity
        intMask = roi > (params.intensityThreshFraction * max(roi(:)));

        % Only edges + intensity
        BW = edges | intMask;
        
        % --- Morphological cleanup ---
        % Remove tiny fragments, close small gaps, optionally fill holes
        BW = bwareaopen(BW, params.shapeMinArea);                      % remove small objects
        BW = imclose(BW, strel('disk', max(1, params.shapeCloseRadius))); % smooth edges
        if params.shapeFillHoles
            BW = imfill(BW,'holes');                                   % fill internal holes
        end

        % --- Label and measure region properties ---
        labels = bwlabel(BW);  % assign unique labels to each connected shape
        props  = regionprops(labels, 'Area','Eccentricity','Solidity','Centroid','Perimeter');
        for n = 1:numel(props)
            % Classify shape based on morphology
            if props(n).Eccentricity < 0.61
                props(n).Class = "circular/elliptical";
            elseif props(n).Solidity < 0.7
                props(n).Class = "Y-shaped/branched";
            else
                props(n).Class = "filament-like";
            end
        end

        % --- Extract global boundaries ---
        boundaries = bwboundaries(BW);
        boundariesGlobal = cell(size(boundaries));
        for n = 1:numel(boundaries)
            b = boundaries{n};
            % Convert from local crop coordinates to full image coordinates
            b(:,1) = b(:,1) + yStart - 1;
            b(:,2) = b(:,2) + xStart - 1;
            boundariesGlobal{n} = b;
        end

        % --- Store results ---
        results(k).BW         = BW;
        results(k).labels     = labels;
        results(k).props      = props;
        results(k).boundaries = boundariesGlobal;
    end
end
