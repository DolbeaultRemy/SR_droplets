function [bottleneckResults] = comprehensiveBottleneckAnalysis(BW, patchProps, thresholdRatio)
% 综合瓶颈分析
% thresholdRatio: 判断为瓶颈的阈值比例（如0.3表示宽度小于最大宽度的30%）

    if nargin < 3
        thresholdRatio = 0.3;
    end
    
    numRegions = length(patchProps);
    bottleneckResults = struct();
    
    for i = 1:numRegions
        regionBW = false(size(BW));
        regionBW(patchProps(i).PixelIdxList) = true;
        
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
        else
            bottleneckResults(i).RegionIndex = i;
            bottleneckResults(i).MinWidth = 0;
            bottleneckResults(i).MaxWidth = 0;
            bottleneckResults(i).MeanWidth = 0;
            bottleneckResults(i).HasBottleneck = false;
            bottleneckResults(i).BottleneckRatio = 1;
            bottleneckResults(i).Area = stats.Area;
        end
    end
end