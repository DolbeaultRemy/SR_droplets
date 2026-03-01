function [bottleneckResults, BW_split, CC_split, splitInfo] = comprehensiveBottleneckAnalysisWithSplitting(BW, CC, thresholdRatio, minArea, maxIterations)
% 综合瓶颈分析，并在瓶颈处分割区域，递归分割直到没有瓶颈点
% 输入：
%   BW - 二值图像
%   CC - bwconncomp的输出
%   thresholdRatio - 判断为瓶颈的阈值比例（如0.3表示宽度小于最大宽度的30%）
%   minArea - 分割后区域的最小面积（可选，默认50）
%   maxIterations - 最大分割迭代次数（可选，默认5）
% 输出：
%   bottleneckResults - 瓶颈分析结果
%   BW_split - 分割后的二值图像
%   CC_split - 分割后的连通组件
%   splitInfo - 分割信息结构体

    if nargin < 3
        thresholdRatio = 0.3;
    end
    if nargin < 4
        minArea = 50; % 默认最小面积
    end
    if nargin < 5
        maxIterations = 5; % 默认最大迭代次数
    end
    
    % 初始化
    BW_current = BW;
    CC_current = CC;
    
    % 存储分割历史
    splitHistory = struct();
    splitHistory(1).BW = BW_current;
    splitHistory(1).CC = CC_current;
    splitHistory(1).iteration = 1;
    
    iteration = 1;
    hasBottleneck = true;
    
    % fprintf('开始瓶颈分析，最大迭代次数: %d\n', maxIterations);
    
    while hasBottleneck && iteration <= maxIterations
        % fprintf('迭代 %d: 分析 %d 个区域...\n', iteration, CC_current.NumObjects);
        
        % 分析当前迭代的瓶颈
        [currentResults, BW_new, CC_new, currentSplitInfo] = ...
            singleIterationBottleneckAnalysis(BW_current, CC_current, thresholdRatio, minArea);
        
        % 检查是否还有瓶颈区域
        hasBottleneck = any([currentResults.HasBottleneck] & [currentResults.Area] > minArea * 2);
        
        % 更新当前状态
        BW_current = BW_new;
        CC_current = CC_new;
        
        % 记录分割历史
        iteration = iteration + 1;
        splitHistory(iteration).BW = BW_current;
        splitHistory(iteration).CC = CC_current;
        splitHistory(iteration).iteration = iteration;
        splitHistory(iteration).results = currentResults;
        splitHistory(iteration).splitInfo = currentSplitInfo;
        
        % fprintf('迭代 %d 完成: 分割出 %d 个新区域，仍有瓶颈: %d\n', ...
        %     iteration-1, currentSplitInfo.numSplits, hasBottleneck);
        
        % 如果没有区域变化，提前退出
        if CC_current.NumObjects == splitHistory(iteration-1).CC.NumObjects
            % fprintf('区域数量未变化，提前退出迭代\n');
            break;
        end
    end
    
    % 最终结果
    BW_split = BW_current;
    CC_split = CC_current;
    
    % 生成最终瓶颈分析结果
    bottleneckResults = analyzeFinalBottlenecks(BW_split, CC_split, thresholdRatio);
    
    % 整理分割信息
    splitInfo = struct();
    splitInfo.splitHistory = splitHistory;
    splitInfo.totalIterations = iteration - 1;
    splitInfo.finalRegionCount = CC_split.NumObjects;
    splitInfo.originalRegionCount = CC.NumObjects;
    splitInfo.maxIterations = maxIterations;
end

function [bottleneckResults, BW_split, CC_split, splitInfo] = singleIterationBottleneckAnalysis(BW, CC, thresholdRatio, minArea)
% 单次迭代的瓶颈分析
    
    numRegions = CC.NumObjects;
    bottleneckResults = struct();
    BW_split = BW; % 初始化为原始图像
    splitInfo = struct();
    
    % 初始化分割信息
    splitInfo.splitPoints = cell(numRegions, 1);
    splitInfo.numSplits = 0;
    splitInfo.regionsSplit = [];
    
    for i = 1:numRegions
        regionBW = false(size(BW));
        regionBW(CC.PixelIdxList{i}) = true;
        
        % 计算区域属性
        stats = regionprops(regionBW, 'BoundingBox', 'Area', 'Perimeter');
        
        % 计算距离变换
        distTransform = bwdist(~regionBW);
        
        % 计算骨架
        skeleton = bwmorph(regionBW, 'thin', inf);
        
        % 获取骨架上的宽度
        skeletonWidths = 2 * distTransform(skeleton);
        
        if ~isempty(skeletonWidths) && stats.Area > minArea
            minWidth = min(skeletonWidths);
            maxWidth = max(skeletonWidths);
            meanWidth = mean(skeletonWidths);
            
            % 判断是否有瓶颈
            isBottleneck = (minWidth / maxWidth) < thresholdRatio;
            
            % 存储结果
            bottleneckResults(i).RegionIndex = i;
            bottleneckResults(i).MinWidth = minWidth;
            bottleneckResults(i).MaxWidth = maxWidth;
            bottleneckResults(i).MeanWidth = meanWidth;
            bottleneckResults(i).HasBottleneck = isBottleneck;
            bottleneckResults(i).BottleneckRatio = minWidth / maxWidth;
            bottleneckResults(i).Area = stats.Area;
            bottleneckResults(i).Skeleton = skeleton;
            bottleneckResults(i).DistTransform = distTransform;
            
            % 如果有瓶颈，进行分割
            if isBottleneck && stats.Area > minArea * 2
                [regionBW_split, splitPoints] = splitRegionAtBottleneck(regionBW, skeleton, distTransform, thresholdRatio);
                
                % 检查分割后的区域数量
                CC_temp = bwconncomp(regionBW_split);
                if CC_temp.NumObjects > 1 % 确保确实分割了
                    % 更新分割后的图像
                    BW_split(CC.PixelIdxList{i}) = false; % 移除原区域
                    BW_split = BW_split | regionBW_split; % 添加分割后的区域
                    
                    % 记录分割信息
                    splitInfo.splitPoints{i} = splitPoints;
                    splitInfo.numSplits = splitInfo.numSplits + 1;
                    splitInfo.regionsSplit(end+1) = i;
                    
                    bottleneckResults(i).WasSplit = true;
                    bottleneckResults(i).NumSplitParts = CC_temp.NumObjects;
                else
                    bottleneckResults(i).WasSplit = false;
                    bottleneckResults(i).NumSplitParts = 1;
                    splitInfo.splitPoints{i} = [];
                end
            else
                bottleneckResults(i).WasSplit = false;
                bottleneckResults(i).NumSplitParts = 1;
                splitInfo.splitPoints{i} = [];
            end
        else
            % 区域太小或没有骨架点
            bottleneckResults(i).RegionIndex = i;
            bottleneckResults(i).MinWidth = 0;
            bottleneckResults(i).MaxWidth = 0;
            bottleneckResults(i).MeanWidth = 0;
            bottleneckResults(i).HasBottleneck = false;
            bottleneckResults(i).BottleneckRatio = 1;
            bottleneckResults(i).Area = stats.Area;
            bottleneckResults(i).WasSplit = false;
            bottleneckResults(i).NumSplitParts = 1;
            splitInfo.splitPoints{i} = [];
        end
    end
    
    % 更新分割后的连通组件
    CC_split = bwconncomp(BW_split);
end

function bottleneckResults = analyzeFinalBottlenecks(BW, CC, thresholdRatio)
% 分析最终图像的瓶颈情况
    
    numRegions = CC.NumObjects;
    bottleneckResults = struct();
    
    for i = 1:numRegions
        regionBW = false(size(BW));
        regionBW(CC.PixelIdxList{i}) = true;
        
        % 计算区域属性
        stats = regionprops(regionBW, 'BoundingBox', 'Area', 'Perimeter');
        
        % 计算距离变换
        distTransform = bwdist(~regionBW);
        
        % 计算骨架
        skeleton = bwmorph(regionBW, 'thin', inf);
        
        % 获取骨架上的宽度
        skeletonWidths = 2 * distTransform(skeleton);
        
        if ~isempty(skeletonWidths)
            minWidth = min(skeletonWidths);
            maxWidth = max(skeletonWidths);
            meanWidth = mean(skeletonWidths);
            
            % 判断是否有瓶颈
            isBottleneck = (minWidth / maxWidth) < thresholdRatio;
            
            % 存储结果
            bottleneckResults(i).RegionIndex = i;
            bottleneckResults(i).MinWidth = minWidth;
            bottleneckResults(i).MaxWidth = maxWidth;
            bottleneckResults(i).MeanWidth = meanWidth;
            bottleneckResults(i).HasBottleneck = isBottleneck;
            bottleneckResults(i).BottleneckRatio = minWidth / maxWidth;
            bottleneckResults(i).Area = stats.Area;
            bottleneckResults(i).Skeleton = skeleton;
            bottleneckResults(i).DistTransform = distTransform;
            bottleneckResults(i).WasSplit = false; % 最终结果不记录分割历史
            bottleneckResults(i).NumSplitParts = 1;
        else
            bottleneckResults(i).RegionIndex = i;
            bottleneckResults(i).MinWidth = 0;
            bottleneckResults(i).MaxWidth = 0;
            bottleneckResults(i).MeanWidth = 0;
            bottleneckResults(i).HasBottleneck = false;
            bottleneckResults(i).BottleneckRatio = 1;
            bottleneckResults(i).Area = stats.Area;
            bottleneckResults(i).WasSplit = false;
            bottleneckResults(i).NumSplitParts = 1;
        end
    end
end

function [regionBW_split, splitPoints] = splitRegionAtBottleneck(regionBW, skeleton, distTransform, thresholdRatio)
% 在瓶颈处分割区域
% 输入：
%   regionBW - 区域二值图像
%   skeleton - 区域骨架
%   distTransform - 区域距离变换
%   minWidth - 瓶颈最小宽度
% 输出：
%   regionBW_split - 分割后的区域二值图像
%   splitPoints - 分割点坐标

    % 找到瓶颈点（最小宽度的骨架点）
    [skeletonY, skeletonX] = find(skeleton);
    widthsAtSkeleton = 2 * distTransform(skeleton);
    maxWidth = max(widthsAtSkeleton);

    flagSplit = false;
    for i=1:10

        minWidth = min(widthsAtSkeleton);
        if (minWidth / maxWidth) >= thresholdRatio
            regionBW_split = regionBW;
            splitPoints = [];
            return
        end
        bottleneckIndices = find(abs(widthsAtSkeleton - minWidth) < 0.1); % 容差
        
        if isempty(bottleneckIndices)
            regionBW_split = regionBW;
            splitPoints = [];
            return;
        end
        
        % 选择最中心的瓶颈点作为分割点
        [centerY, centerX] = find(regionBW);
        centerY = mean(centerY);
        centerX = mean(centerX);
        
        bottleneckX = skeletonX(bottleneckIndices);
        bottleneckY = skeletonY(bottleneckIndices);
        
        % 计算每个瓶颈点到区域中心的距离
        distancesToCenter = sqrt((bottleneckX - centerX).^2 + (bottleneckY - centerY).^2);
        [~, closestIdx] = min(distancesToCenter);
        
        splitX = bottleneckX(closestIdx);
        splitY = bottleneckY(closestIdx);
        splitPoints = [splitX, splitY];
        
        % 方法1: 使用形态学操作在瓶颈处断开
        regionBW_split = splitWithMorphology(regionBW, splitX, splitY, minWidth);

        [~, num] = bwlabel(regionBW_split); 

        if num>1
            return;
        else
            widthsAtSkeleton(bottleneckIndices) = maxWidth;
        end
    
    end
    
    % 方法2: 如果方法1无效，使用更激进的分割方法
    % if length(unique(regionBW_split(regionBW_split))) <= 2 % 如果没有成功分割
    %     regionBW_split = splitWithLineCut(regionBW, splitX, splitY, minWidth);
    % end
end

function regionBW_split = splitWithMorphology(regionBW, splitX, splitY, minWidth)
% 使用形态学操作在指定点分割区域
    
    % 创建分割结构元素（大小基于瓶颈宽度）
    seSize = max(1, ceil(double(minWidth) / 2));
    se = strel('disk', seSize);
    
    % 创建分割掩码
    splitMask = false(size(regionBW));
    splitMask(splitY, splitX) = true;
    
    % 膨胀分割点
    splitMask = imdilate(splitMask, se);
    
    % 从原区域中移除分割点区域
    regionBW_split = regionBW & ~splitMask;
    
    % 清理小的孤立区域
    regionBW_split = bwareaopen(regionBW_split, ceil(double(minWidth)^2));
end

function originalToNewMap = createRegionMapping(BW_original, BW_split, CC_original, CC_split)
% 建立原始区域到新区域的映射（简化版）
    originalToNewMap = cell(CC_original.NumObjects, 1);
    
    for i = 1:CC_original.NumObjects
        originalPixels = CC_original.PixelIdxList{i};
        overlappingRegions = [];
        
        for j = 1:CC_split.NumObjects
            splitPixels = CC_split.PixelIdxList{j};
            % 检查是否有重叠
            if any(ismember(originalPixels, splitPixels))
                overlappingRegions(end+1) = j;
            end
        end
        
        originalToNewMap{i} = overlappingRegions;
    end
end