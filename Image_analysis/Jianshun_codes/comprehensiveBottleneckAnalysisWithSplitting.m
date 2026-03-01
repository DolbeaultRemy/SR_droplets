function [bottleneckResults, BW_split, CC_split, splitInfo] = ...
    comprehensiveBottleneckAnalysisWithSplitting(BW, CC, thresholdRatio, minArea, maxIterations)
% COMPREHENSIVEBOTTLENECKANALYSISWITHSPLITTING
%   Recursively split regions at bottlenecks and analyze final regions.
% Inputs:
%   BW             - binary image
%   CC             - bwconncomp output of BW
%   thresholdRatio - width ratio below which a bottleneck is considered (e.g., 0.3)
%   minArea        - minimum area to consider splitting (optional, default 50)
%   maxIterations  - maximum number of splitting iterations (optional, default 5)
% Outputs:
%   bottleneckResults - struct array with analysis for each final region
%   BW_split          - binary image after splitting
%   CC_split          - bwconncomp of BW_split
%   splitInfo         - struct with splitting history

    if nargin < 3 || isempty(thresholdRatio); thresholdRatio = 0.3; end
    if nargin < 4 || isempty(minArea); minArea = 50; end
    if nargin < 5 || isempty(maxIterations); maxIterations = 5; end

    BW_current = BW;
    CC_current = CC;
    splitHistory = struct();
    splitHistory(1).BW = BW_current;
    splitHistory(1).CC = CC_current;
    splitHistory(1).iteration = 1;

    iteration = 1;
    hasBottleneck = true;

    while hasBottleneck && iteration <= maxIterations
        [currentResults, BW_new, CC_new, currentSplitInfo] = ...
            singleIterationBottleneckAnalysis(BW_current, CC_current, thresholdRatio, minArea);

        % Check if any region still has a bottleneck (and area large enough)
        hasBottleneck = any([currentResults.HasBottleneck] & [currentResults.Area] > minArea * 2);

        BW_current = BW_new;
        CC_current = CC_new;

        iteration = iteration + 1;
        splitHistory(iteration).BW = BW_current;
        splitHistory(iteration).CC = CC_current;
        splitHistory(iteration).iteration = iteration;
        splitHistory(iteration).results = currentResults;
        splitHistory(iteration).splitInfo = currentSplitInfo;

        % Stop if region count did not change
        if CC_current.NumObjects == splitHistory(iteration-1).CC.NumObjects
            break;
        end
    end

    BW_split = BW_current;
    CC_split = CC_current;
    bottleneckResults = analyzeFinalBottlenecks(BW_split, CC_split, thresholdRatio);

    splitInfo = struct();
    splitInfo.splitHistory = splitHistory;
    splitInfo.totalIterations = iteration - 1;
    splitInfo.finalRegionCount = CC_split.NumObjects;
    splitInfo.originalRegionCount = CC.NumObjects;
    splitInfo.maxIterations = maxIterations;
end

%% ------------------------------------------------------------------------
function [bottleneckResults, BW_split, CC_split, splitInfo] = ...
    singleIterationBottleneckAnalysis(BW, CC, thresholdRatio, minArea)
% Single iteration of bottleneck analysis: identify bottlenecks and split.
    numRegions = CC.NumObjects;
    bottleneckResults = struct();
    BW_split = BW; % start with original
    splitInfo = struct('splitPoints', cell(numRegions,1), 'numSplits', 0, 'regionsSplit', []);

    for i = 1:numRegions
        regionBW = false(size(BW));
        regionBW(CC.PixelIdxList{i}) = true;
        stats = regionprops(regionBW, 'BoundingBox', 'Area', 'Perimeter');

        % Distance transform and skeleton
        distTransform = bwdist(~regionBW);
        skeleton = bwmorph(regionBW, 'thin', inf);
        skeletonWidths = 2 * distTransform(skeleton);

        if ~isempty(skeletonWidths) && stats.Area > minArea
            minW = min(skeletonWidths);
            maxW = max(skeletonWidths);
            meanW = mean(skeletonWidths);
            isBottleneck = (minW / maxW) < thresholdRatio;

            % Store results
            res = struct('RegionIndex', i, 'MinWidth', minW, 'MaxWidth', maxW, ...
                'MeanWidth', meanW, 'HasBottleneck', isBottleneck, ...
                'BottleneckRatio', minW/maxW, 'Area', stats.Area, ...
                'Skeleton', skeleton, 'DistTransform', distTransform);

            if isBottleneck && stats.Area > minArea * 2
                [regionBW_split, splitPoints] = splitRegionAtBottleneck(regionBW, skeleton, distTransform, thresholdRatio);
                CC_temp = bwconncomp(regionBW_split);
                if CC_temp.NumObjects > 1
                    BW_split(CC.PixelIdxList{i}) = false;
                    BW_split = BW_split | regionBW_split;
                    splitInfo.splitPoints{i} = splitPoints;
                    splitInfo.numSplits = splitInfo.numSplits + 1;
                    splitInfo.regionsSplit(end+1) = i;
                    res.WasSplit = true;
                    res.NumSplitParts = CC_temp.NumObjects;
                else
                    res.WasSplit = false;
                    res.NumSplitParts = 1;
                end
            else
                res.WasSplit = false;
                res.NumSplitParts = 1;
            end
            bottleneckResults(i) = res;
        else
            % Region too small or no skeleton
            bottleneckResults(i) = struct('RegionIndex', i, 'MinWidth', 0, 'MaxWidth', 0, ...
                'MeanWidth', 0, 'HasBottleneck', false, 'BottleneckRatio', 1, ...
                'Area', stats.Area, 'WasSplit', false, 'NumSplitParts', 1);
        end
    end
    CC_split = bwconncomp(BW_split);
end

%% ------------------------------------------------------------------------
function bottleneckResults = analyzeFinalBottlenecks(BW, CC, thresholdRatio)
% Analyze final regions after splitting.
    numRegions = CC.NumObjects;
    bottleneckResults = struct();
    for i = 1:numRegions
        regionBW = false(size(BW));
        regionBW(CC.PixelIdxList{i}) = true;
        stats = regionprops(regionBW, 'Area');
        distTransform = bwdist(~regionBW);
        skeleton = bwmorph(regionBW, 'thin', inf);
        skeletonWidths = 2 * distTransform(skeleton);

        if ~isempty(skeletonWidths)
            minW = min(skeletonWidths);
            maxW = max(skeletonWidths);
            meanW = mean(skeletonWidths);
            isBottleneck = (minW / maxW) < thresholdRatio;
            bottleneckResults(i) = struct('RegionIndex', i, 'MinWidth', minW, ...
                'MaxWidth', maxW, 'MeanWidth', meanW, 'HasBottleneck', isBottleneck, ...
                'BottleneckRatio', minW/maxW, 'Area', stats.Area, ...
                'WasSplit', false, 'NumSplitParts', 1);
        else
            bottleneckResults(i) = struct('RegionIndex', i, 'MinWidth', 0, ...
                'MaxWidth', 0, 'MeanWidth', 0, 'HasBottleneck', false, ...
                'BottleneckRatio', 1, 'Area', stats.Area, ...
                'WasSplit', false, 'NumSplitParts', 1);
        end
    end
end

%% ------------------------------------------------------------------------
function [regionBW_split, splitPoints] = splitRegionAtBottleneck(regionBW, skeleton, distTransform, thresholdRatio)
% Attempt to split region at bottleneck point(s).
    [skeletonY, skeletonX] = find(skeleton);
    widthsAtSkeleton = 2 * distTransform(skeleton);
    maxWidth = max(widthsAtSkeleton);

    % Iteratively try to split at the narrowest point until success or limit
    for attempt = 1:10
        minWidth = min(widthsAtSkeleton);
        if (minWidth / maxWidth) >= thresholdRatio
            regionBW_split = regionBW;
            splitPoints = [];
            return;
        end
        bottleneckIndices = find(abs(widthsAtSkeleton - minWidth) < 0.1);
        if isempty(bottleneckIndices)
            regionBW_split = regionBW;
            splitPoints = [];
            return;
        end

        % Choose bottleneck point closest to region center
        [centerY, centerX] = find(regionBW);
        centerY = mean(centerY);
        centerX = mean(centerX);
        bottleneckX = skeletonX(bottleneckIndices);
        bottleneckY = skeletonY(bottleneckIndices);
        distToCenter = sqrt((bottleneckX - centerX).^2 + (bottleneckY - centerY).^2);
        [~, closestIdx] = min(distToCenter);
        splitX = bottleneckX(closestIdx);
        splitY = bottleneckY(closestIdx);
        splitPoints = [splitX, splitY];

        % Morphological split: remove a disk at the bottleneck
        regionBW_split = splitWithMorphology(regionBW, splitX, splitY, minWidth);

        [~, num] = bwlabel(regionBW_split);
        if num > 1
            return; % success
        else
            % Mark this bottleneck as "not narrow" and try next narrowest
            widthsAtSkeleton(bottleneckIndices) = maxWidth;
        end
    end
    % If all attempts fail, return original
    regionBW_split = regionBW;
    splitPoints = [];
end

%% ------------------------------------------------------------------------
function regionBW_split = splitWithMorphology(regionBW, splitX, splitY, minWidth)
% Split region by removing a disk of radius based on minWidth.
    seSize = max(1, ceil(double(minWidth) / 2));
    se = strel('disk', seSize);
    splitMask = false(size(regionBW));
    splitMask(splitY, splitX) = true;
    splitMask = imdilate(splitMask, se);
    regionBW_split = regionBW & ~splitMask;
    regionBW_split = bwareaopen(regionBW_split, ceil(double(minWidth)^2));
end