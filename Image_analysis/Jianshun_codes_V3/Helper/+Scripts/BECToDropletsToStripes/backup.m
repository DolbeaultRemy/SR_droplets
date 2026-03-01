% function [bottleneckResults, BW_split, CC_split, splitInfo] = comprehensiveBottleneckAnalysisWithSplitting(BW, CC, thresholdRatio, minArea)
% % 综合瓶颈分析，并在瓶颈处分割区域
% % 输入：
% %   BW - 二值图像
% %   CC - bwconncomp的输出
% %   thresholdRatio - 判断为瓶颈的阈值比例（如0.3表示宽度小于最大宽度的30%）
% %   minArea - 分割后区域的最小面积（可选，默认50）
% % 输出：
% %   bottleneckResults - 瓶颈分析结果
% %   BW_split - 分割后的二值图像
% %   CC_split - 分割后的连通组件
% %   splitInfo - 分割信息结构体
% 
%     if nargin < 3
%         thresholdRatio = 0.3;
%     end
%     if nargin < 4
%         minArea = 50; % 默认最小面积
%     end
% 
%     numRegions = CC.NumObjects;
%     bottleneckResults = struct();
%     BW_split = BW; % 初始化为原始图像
%     splitInfo = struct();
% 
%     % 初始化分割信息
%     splitInfo.originalToNewMap = cell(numRegions, 1);
%     splitInfo.splitPoints = cell(numRegions, 1);
%     splitInfo.numSplits = 0;
% 
%     for i = 1:numRegions
%         regionBW = false(size(BW));
%         regionBW(CC.PixelIdxList{i}) = true;
% 
%         % 计算区域属性
%         stats = regionprops(regionBW, 'BoundingBox', 'Area', 'Perimeter');
% 
%         % 计算距离变换
%         distTransform = bwdist(~regionBW);
% 
%         % 计算骨架
%         skeleton = bwmorph(regionBW, 'thin', inf);
% 
%         % 获取骨架上的宽度
%         skeletonWidths = 2 * distTransform(skeleton);
% 
%         if ~isempty(skeletonWidths)
%             minWidth = min(skeletonWidths);
%             maxWidth = max(skeletonWidths);
%             meanWidth = mean(skeletonWidths);
% 
%             % 判断是否有瓶颈
%             isBottleneck = (minWidth / maxWidth) < thresholdRatio;
% 
%             % 存储结果
%             bottleneckResults(i).RegionIndex = i;
%             bottleneckResults(i).MinWidth = minWidth;
%             bottleneckResults(i).MaxWidth = maxWidth;
%             bottleneckResults(i).MeanWidth = meanWidth;
%             bottleneckResults(i).HasBottleneck = isBottleneck;
%             bottleneckResults(i).BottleneckRatio = minWidth / maxWidth;
%             bottleneckResults(i).Area = stats.Area;
%             bottleneckResults(i).Skeleton = skeleton;
%             bottleneckResults(i).DistTransform = distTransform;
% 
%             % 如果有瓶颈，进行分割
%             if isBottleneck && stats.Area > minArea * 2 % 只有面积足够大才分割
%                 [regionBW_split, splitPoints] = splitRegionAtBottleneck(regionBW, skeleton, distTransform, thresholdRatio);
% 
%                 % 更新分割后的图像
%                 BW_split(CC.PixelIdxList{i}) = false; % 移除原区域
%                 BW_split = BW_split | regionBW_split; % 添加分割后的区域
% 
%                 % 记录分割信息
%                 splitInfo.splitPoints{i} = splitPoints;
%                 splitInfo.numSplits = splitInfo.numSplits + 1;
% 
%                 bottleneckResults(i).WasSplit = true;
%                 bottleneckResults(i).NumSplitParts = length(unique(regionBW_split(regionBW_split))) - 1; % 减去背景
%             else
%                 bottleneckResults(i).WasSplit = false;
%                 bottleneckResults(i).NumSplitParts = 1;
%                 splitInfo.splitPoints{i} = [];
%             end
%         else
%             bottleneckResults(i).RegionIndex = i;
%             bottleneckResults(i).MinWidth = 0;
%             bottleneckResults(i).MaxWidth = 0;
%             bottleneckResults(i).MeanWidth = 0;
%             bottleneckResults(i).HasBottleneck = false;
%             bottleneckResults(i).BottleneckRatio = 1;
%             bottleneckResults(i).Area = stats.Area;
%             bottleneckResults(i).WasSplit = false;
%             bottleneckResults(i).NumSplitParts = 1;
%             splitInfo.splitPoints{i} = [];
%         end
%     end
% 
%     % 更新分割后的连通组件
%     CC_split = bwconncomp(BW_split);
% 
%     % 建立原始区域到新区域的映射
%     splitInfo.originalToNewMap = createRegionMapping(BW, BW_split, CC, CC_split);
% end

% function originalToNewMap = createRegionMapping(originalBW, splitBW, originalCC, splitCC)
% % 建立原始区域到分割后区域的映射
% 
%     numOriginal = originalCC.NumObjects;
%     originalToNewMap = cell(numOriginal, 1);
% 
%     % 创建原始标签矩阵
%     originalLabels = labelmatrix(originalCC);
% 
%     % 创建分割后标签矩阵
%     splitLabels = labelmatrix(splitCC);
% 
%     for i = 1:numOriginal
%         % 获取原始区域像素
%         originalPixels = originalCC.PixelIdxList{i};
% 
%         % 找出这些像素在分割后图像中属于哪些区域
%         splitRegions = splitLabels(originalPixels);
%         splitRegions = splitRegions(splitRegions > 0); % 去除背景
% 
%         if isempty(splitRegions)
%             originalToNewMap{i} = [];
%         else
%             % 找出所有相关的分割区域
%             uniqueRegions = unique(splitRegions);
% 
%             % 计算每个分割区域与原始区域的重叠比例
%             overlapRatios = zeros(size(uniqueRegions));
%             for j = 1:length(uniqueRegions)
%                 regionPixels = splitCC.PixelIdxList{uniqueRegions(j)};
%                 overlap = sum(ismember(regionPixels, originalPixels));
%                 overlapRatios(j) = overlap / length(regionPixels);
%             end
% 
%             % 只保留重叠比例足够高的区域
%             validRegions = uniqueRegions(overlapRatios > 0.5);
%             originalToNewMap{i} = validRegions;
%         end
%     end
% end