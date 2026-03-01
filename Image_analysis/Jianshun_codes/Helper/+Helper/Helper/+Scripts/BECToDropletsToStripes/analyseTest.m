% folderPath = "C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2025\10\14\0003\AnalysisPlot\Data\";
% folderPath = "C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2026\01\01\0005\AnalysisPlot\Data\";
folderPath = "D:\Jianshun\Data\2026\01\12\0002\AnalysisPlot\Data\";
file_list = dir(fullfile(folderPath, '*.mat'));

load(fullfile(folderPath, file_list(1).name));

plotIdx = length(file_list);
plotFlag = false;

number = zeros(length(scan_parameter_values), 1);
averageArea = zeros(length(scan_parameter_values), 1);
averageSumOD = zeros(length(scan_parameter_values), 1);
sumOD = zeros(length(scan_parameter_values), 150, 1);
detectedArea = zeros(length(scan_parameter_values), 150, 1);
averageOD = zeros(length(scan_parameter_values), 150, 1);
selectFlag = zeros(length(scan_parameter_values), 150, 1);
aspectRatio = zeros(length(scan_parameter_values), 150, 1);
averageAspectRatio = zeros(length(scan_parameter_values), 1);
allRotatedPoints = zeros(length(scan_parameter_values), 150, 2);
edgeLengths = zeros(length(scan_parameter_values), 150, 1);
averageEdgeLengths = zeros(length(scan_parameter_values), 1);
allAngles = zeros(length(scan_parameter_values), 5000);
% allRotatedPointSingle = [];

for i = 1:plotIdx

    load(fullfile(folderPath, file_list(i).name));

    if CC_split.NumObjects == 0
        averageArea(i) = 0;
    else
        % 计算每个区域的面积（像素数量）
        areas = cellfun(@numel, CC_split.PixelIdxList);
        detectedArea(i,1:(length(areas))) = areas;
    end

    numRegions = length(CC_split.PixelIdxList);

    % 遍历每个区域，计算灰度值总和
    for k = 1:numRegions
        % 获取当前区域的像素索引
        pixelIdx = cell2mat(CC_split.PixelIdxList(k));
    
        % 计算该区域的灰度值总和
        sumOD(i, k) = sum(imgCropped(pixelIdx));

        averageOD(i, k) = sumOD(i, k) / detectedArea(i, k);
    end

    for j = 1:CC_split.NumObjects
        if averageOD(i, j) > 2.5 && detectedArea(i, j) > 30
            selectFlag(i, j) = 1;
        end
    end

    aspectRatioTemp = calculate_aspect_ratios_simple(CC_split);
    aspectRatio(i, 1:length(aspectRatioTemp)) = aspectRatioTemp;

    for j = 1:CC_split.NumObjects
        if selectFlag(i, j)
            averageArea(i) = averageArea(i) + detectedArea(i, j);
            averageSumOD(i) = averageSumOD(i) + sumOD(i, j);
            averageAspectRatio(i) = averageAspectRatio(i) + aspectRatio(i, j);
        end
    end

    number(i) = sum(selectFlag(i, :));
    averageArea(i) = averageArea(i) / number(i);
    averageSumOD(i) = averageSumOD(i) / number(i);
    averageAspectRatio(i) = averageAspectRatio(i) / number(i);

    centers = calculateWeightedCenters(CC_split, imgCropped, selectFlag(i, :));
    
    [allRotatedPointSingle, edgeLength, edgesPlot] = simpleAlignedNewtonDiagram(centers);
    if length(allRotatedPointSingle)>0
        allRotatedPoints(i, 1:length(allRotatedPointSingle(:, 1)), :) = allRotatedPointSingle;
    end
    edgeLengths(i, 1:length(edgeLength), :) = edgeLength;
    averageEdgeLengths(i) = mean(edgeLength);

    [allAngle, edgeLengths, edgesPlot] = edgeAngleHistogram(centers);
    if length(allAngle)>0
        allAngles(i, 1:length(allAngle(:, 1))) = allAngle;
    end
    
    % for j = 1:CC_split.NumObjects
    %     if selectFlag(i, j)
    %         regionBW = false(size(imgCropped));
    %         regionBW(CC.PixelIdxList(j)) = true;
    % 
    %         % 计算距离变换
    %         distTransform = bwdist(~regionBW);
    % 
    %         % 计算骨架
    %         skeleton = bwmorph(regionBW, 'thin', inf);
    % 
    %         % 获取骨架上的宽度
    %         skeletonWidths = 2 * distTransform(skeleton);
    %     end
    % end
    
    if plotFlag
        mask = false(size(imgCropped));
    
        % 将CC_split中的所有区域合并到mask中
        for k = 1:CC_split.NumObjects
            if selectFlag(i, k)
                mask(CC_split.PixelIdxList{k}) = true;
            end
        end
    
        % 获取区域的边界
        boundaries = bwboundaries(mask);
    
        figure(2)
    
        % 显示原始图像
        imagesc(imgCropped, [0, 3])
        colormap jet
        clim([-0.5 3.5])
        hold on;
    
        % 绘制所有区域的边界
        for k = 1:length(boundaries)
            boundary = boundaries{k};
            plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2);
        end
    
        % 创建与原始图像相同大小的空白蒙版
        try
            DT = delaunayTriangulation(centers);
    
            % 获取凸包上的点索引（最外层点）
            convexHullIndices = convexHull(DT);
        catch
            title("B = " + sprintf('%.2f (G)', scan_parameter_values(i)))
    
            theme("light")
            set(gca, 'FontSize', 16);
            set(gcf,'Position',[-1200 200 900 600])
    
            hold off;
    
            save_path = folderPath + "FinialPlot\";
            if ~exist(save_path, 'dir')
                mkdir(save_path);
            end
    
            filename = int2str(int64(i)) + "_B = " + sprintf('%.2f (G)', scan_parameter_values(i)) + ".jpg";
            full_path = fullfile(save_path, filename);
            % 导出图像
            print(gcf, full_path, '-djpeg', '-r300');
    
            continue
        end
    
        % 注意：convexHull 返回的是凸包多边形的顶点索引（首尾相连）
        % 移除重复点（首尾相同）
        if convexHullIndices(1) == convexHullIndices(end)
            convexHullIndices = convexHullIndices(1:end-1);
        end
    
        % 创建逻辑索引标记哪些点是最外层点
        isOuterPoint = false(size(centers, 1), 1);
        isOuterPoint(convexHullIndices) = true;
    
        % 移除最外层点，保留内部点
        innerPoints = centers(~isOuterPoint, :);
        outerPoints = centers(isOuterPoint, :);
    
        % 绘制Delaunay三角网格
        triplot(DT, 'b-', 'LineWidth', 2.5);
    
        nEdges = size(edgesPlot, 1);
        x_red = zeros(3*nEdges, 1);
        y_red = zeros(3*nEdges, 1);
    
        for k = 1:nEdges
            idx = 3*(k-1) + 1;
            x_red(idx:idx+1) = DT.Points(edgesPlot(k, :), 1);
            y_red(idx:idx+1) = DT.Points(edgesPlot(k, :), 2);
            x_red(idx+2) = NaN;
            y_red(idx+2) = NaN;
        end
    
        % 一次性绘制所有红色边
        plot(x_red, y_red, 'r-', 'LineWidth', 3);
    
        % 标记中心点
        plot(innerPoints(:, 1), innerPoints(:, 2), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
        plot(innerPoints(:, 1), innerPoints(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        plot(outerPoints(:, 1), outerPoints(:, 2), 'k+', 'MarkerSize', 15, 'LineWidth', 2);
        plot(outerPoints(:, 1), outerPoints(:, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    
        title("B = " + sprintf('%.2f (G)', scan_parameter_values(i)))
    
        theme("light")
        set(gca, 'FontSize', 16);
        set(gcf,'Position',[-1200 200 900 600])
    
        hold off;
    
        save_path = folderPath + "FinialPlot\";
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
    
        filename = int2str(int64(i)) + "_B = " + sprintf('%.2f (G)', scan_parameter_values(i)) + ".jpg";
        full_path = fullfile(save_path, filename);
        % 导出图像
        print(gcf, full_path, '-djpeg', '-r300'); 
    end

end

aspectRatio(~selectFlag) = 0;

%%

figure

y_data = number;
x_data = scan_parameter_values;

% 分组计算
unique_x = unique(x_data);
means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);

% 绘制
% figure;
% subplot(1, 3, 1);
subplot(2, 2, 1);
errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'MarkerFaceColor', 'blue');
% xlabel('B (G)');
xlabel('\alpha (deg)');
ylabel('Detected Area Number');
% title('分组统计图');
grid on;
theme("light")
set(gca, 'FontSize', 16);
xlim([min(x_data) max(x_data)])
% set(gcf,'Position',[-1200 200 900 600])
% 
% save_path = folderPath;
% filename = "number.jpg";
% full_path = fullfile(save_path, filename);
% print(gcf, full_path, '-djpeg', '-r300');

%%
y_data = averageArea;
x_data = scan_parameter_values;

% 分组计算
unique_x = unique(x_data);
means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);

% 绘制
% figure;
% subplot(1, 3, 2);
subplot(2, 2, 2);
errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'MarkerFaceColor', 'blue');
% xlabel('B (G)');
xlabel('\alpha (deg)');
ylabel('Detected Area Size');
% title('分组统计图');
grid on;
theme("light")
set(gca, 'FontSize', 16);
xlim([min(x_data) max(x_data)])
% set(gcf,'Position',[-1200 200 900 600])
% 
% save_path = folderPath;
% filename = "averageArea.jpg";
% full_path = fullfile(save_path, filename);
% print(gcf, full_path, '-djpeg', '-r300');

%%
% y_data = averageSumOD;
% x_data = scan_parameter_values;
% 
% % 分组计算
% unique_x = unique(x_data);
% means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
% stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);
% 
% % 绘制
% figure;
% errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
%          'MarkerSize', 8, 'MarkerFaceColor', 'blue');
% % xlabel('B (G)');
% xlabel('\alpha (deg)');
% ylabel('Average Sum OD in Detected Area');
% % title('分组统计图');
% grid on;
% theme("light")
% set(gca, 'FontSize', 16);
% set(gcf,'Position',[-1200 200 900 600])

save_path = folderPath;
filename = "averageSumOD.jpg";
full_path = fullfile(save_path, filename);
print(gcf, full_path, '-djpeg', '-r300');

%%
% % 假设已有数据:
% % averageSumOD - 大小为(n,m)的二维数组
% % scan_parameter_values - 长度为n的一维数组
% 
% % 步骤1: 对行按照scan_parameter_values排序
% % 首先获取唯一的scan_parameter_values和对应的索引
% [unique_scan_values, ~, idx] = unique(scan_parameter_values, 'stable');
% 
% % 计算每个唯一值对应的平均索引位置
% group_indices = accumarray(idx, (1:length(scan_parameter_values))', [], @(x) {x});
% 
% % 创建新的排序索引
% sorted_row_indices = [];
% for i = 1:length(group_indices)
%     sorted_row_indices = [sorted_row_indices; group_indices{i}];
% end
% 
% % 按照排序后的索引重新排列数组
% sorted_averageSumOD = averageOD(sorted_row_indices, :);
% sorted_scan_values = scan_parameter_values(sorted_row_indices);
% 
% % 步骤2: 对每列按照从大到小排序
% % 对每一列单独排序（降序）
% [~, col_indices] = sort(sorted_averageSumOD, 1, 'descend');
% 
% % 由于列排序会影响行数据的一致性，这里我们选择对转置后的矩阵按行排序
% % 这样能保持每行数据的完整性
% sorted_averageSumOD_final = zeros(size(sorted_averageSumOD));
% for col = 1:size(sorted_averageSumOD, 2)
%     sorted_averageSumOD_final(:, col) = sorted_averageSumOD(col_indices(:, col), col);
% end
% 
% % 步骤3: 绘制线图
% figure;
% hold on;
% 
% % 为不同的曲线设置不同的颜色
% colors = lines(size(sorted_averageSumOD_final, 1));
% 
% % 绘制每条曲线
% for i = 1:size(sorted_averageSumOD_final, 1)
%     plot(sorted_averageSumOD_final(i, :), ...
%          'Color', colors(i, :), ...
%          'LineWidth', 1.5, ...
%          'DisplayName', sprintf('B = %.2f', sorted_scan_values(i)));
% end
% 
% % 添加图形标签和标题
% xlabel('Detected Area Index');
% ylabel('Average OD in Detected Area');
% % title('按扫描参数排序的 averageSumOD 曲线');
% legend('show', 'Location', 'best');
% grid on;
% 
% % 美化图形
% set(gca, 'FontSize', 16);
% hold off;
% theme("light")
% set(gcf,'Position',[-1200 200 900 600])

%%
% % 假设已有数据:
% % averageSumOD - 大小为(n,m)的二维数组
% % scan_parameter_values - 长度为n的一维数组
% 
% % 假设已有数据:
% % averageSumOD - 大小为(n,m)的二维数组
% % scan_parameter_values - 长度为n的一维数组
% 
% % 步骤1: 对相同scan_parameter_values的行进行分组并计算平均值和标准差
% [unique_scan_values, ~, idx] = unique(scan_parameter_values, 'sorted');
% 
% % 计算每个唯一值对应的平均值和标准差
% average_data = zeros(length(unique_scan_values), size(averageOD, 2));
% std_data = zeros(length(unique_scan_values), size(averageOD, 2));
% 
% for i = 1:length(unique_scan_values)
%     group_rows = averageOD(idx == i, :);
%     average_data(i, :) = mean(group_rows, 1);
%     std_data(i, :) = std(group_rows, 0, 1);
% end
% 
% % 步骤2: 对列按照从大到小排序（基于平均值）
% % 按列均值降序排序
% [~, col_order] = sort(mean(average_data, 1), 'descend');
% 
% average_data_sorted = average_data(:, col_order);
% std_data_sorted = std_data(:, col_order);
% 
% % 步骤3: 绘制带误差棒的线图
% figure;
% hold on;
% 
% % 使用jet colormap为不同曲线分配渐变色
% num_curves = size(average_data_sorted, 1);
% colors = jet(num_curves);
% 
% % 绘制每条曲线（带误差棒）
% x_values = 1:size(average_data_sorted, 2);
% 
% for i = 1:num_curves
%     % 绘制平均值曲线
%     plot(x_values, average_data_sorted(i, :), ...
%          'Color', colors(i, :), ...
%          'LineWidth', 2, ...
%          'Marker', 'o', ...
%          'MarkerSize', 4, ...
%          'MarkerFaceColor', colors(i, :), ...
%          'DisplayName', sprintf('B = %.2f', unique_scan_values(i)));
% 
%     % 添加误差棒
%     errorbar(x_values, average_data_sorted(i, :), std_data_sorted(i, :), ...
%              'Color', colors(i, :), ...
%              'LineStyle', 'none', ...
%              'LineWidth', 1, ...
%              'CapSize', 5);
% end
% 
% % 添加图形标签和标题
% xlabel('Detected Area Index');
% ylabel('Average OD in Detected Area');
% % title('按扫描参数分组的 averageSumOD 平均值曲线 (±标准差)');
% legend('show', 'Location', 'best');
% grid off;
% 
% % 设置x轴刻度为整数
% set(gca, 'XTick', 1:size(average_data_sorted, 2));
% 
% % 美化图形
% set(gca, 'FontSize', 16);
% hold off;
% theme("light")
% set(gcf,'Position',[-1200 200 900 600])
% 
% % 添加colorbar以显示颜色与扫描参数值的对应关系
% colormap(jet);
% c = colorbar;
% c.Label.String = '扫描参数值';
% c.Ticks = linspace(0, 1, 5); % 根据需要调整刻度数量
% c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ...
%                         linspace(min(unique_scan_values), max(unique_scan_values), 5), ...
%                         'UniformOutput', false);

%%
% 1. 对aspectRatio的行以scan_parameter_values为关键字排序
[sorted_scan, sort_idx] = sort(scan_parameter_values);
sorted_aspectRatio = aspectRatio(sort_idx, :);

% 2. 准备数据
n = size(sorted_aspectRatio, 1);
m = size(sorted_aspectRatio, 2);

% 创建对应的x和y数据
x_data = repelem(sorted_scan, m);  % 每个scan参数值重复m次
y_data = sorted_aspectRatio';    % 展平aspectRatio矩阵
y_data = y_data(:);    % 展平aspectRatio矩阵

% 移除NaN或无穷大的数据点
valid_idx = isfinite(y_data) & isfinite(x_data);
x_data = x_data(valid_idx);
y_data = y_data(valid_idx);

% 3. 设置y轴的bin数量（可调参数）
y_bins = 50;  % 你可以调整这个值

% 4. 计算二维直方图（x轴不划分bin，只对y轴划分bin）
% 获取唯一的x值
unique_x = unique(sorted_scan);
n_x = length(unique_x);

% 初始化计数矩阵
counts = zeros(n_x, y_bins);

% 计算y轴的bin边界
y_min = min(y_data);
y_max = max(y_data);
y_edges = linspace(y_min, y_max, y_bins + 1);
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

% 对每个唯一的x值，统计y值的分布
for i = 1:n_x
    current_x = unique_x(i);
    % 找到当前x对应的所有y值
    y_values = y_data(x_data == current_x);
    
    if ~isempty(y_values)
        % 计算直方图
        hist_counts = histcounts(y_values, y_edges);
        counts(i, :) = hist_counts;
    end
end

% 5. 将计数为0的位置设为NaN（不显示）
counts(counts == 0) = NaN;
counts(:, 1) = NaN;

% 6. 创建网格用于绘图
[X, Y] = meshgrid(unique_x, y_centers);

% 7. 绘制图像
% figure;
% subplot(1, 3, 3);
subplot(2, 2, 3);
hold on
h = pcolor(X, Y, counts');
set(h, 'EdgeColor', 'none');

% 8. 设置colormap
colormap(jet);
colorbar;

% 9. 添加标签和标题
% xlabel('B (G)');
xlabel('\alpha (deg)');
ylabel('Aspect Ratio');
title('Aspect Ratio Distribution vs Scan Parameter');

% % 10. 美化图形
% set(gca, 'FontSize', 16);
% grid off;
% theme("light")
% set(gcf,'Position',[-1200 200 900 600])

% 可选：如果数据范围很大，可以使用对数色标
% if max(counts(:)) / min(counts(counts > 0)) > 100
%     set(gca, 'ColorScale', 'log');
%     colorbar;
%     title('Aspect Ratio Distribution vs Scan Parameter (Log Color Scale)');
% end

y_data = averageAspectRatio;
x_data = scan_parameter_values;

% 分组计算
unique_x = unique(x_data);
means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);

% 绘制
% figure;
a = errorbar(unique_x, means, stds, 'ko-', 'LineWidth', 3, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'black', "DisplayName", "Average value");
% xlabel('B (G)');
xlabel('\alpha (deg)');
ylabel('Average Aspect Ratio');
% title('分组统计图');
grid on;
theme("light")
set(gca, 'FontSize', 16);
% set(gcf,'Position',[-1200 200 900 600])
xlim([min(x_data) max(x_data)])
legend(a);
hold off

% set(gcf,'Position',[-1500 200 600*2 400*2])

save_path = folderPath;
filename = "Result.jpg";
full_path = fullfile(save_path, filename);
print(gcf, full_path, '-djpeg', '-r300');

% save_path = folderPath;
% filename = "aspectRatio.jpg";
% full_path = fullfile(save_path, filename);
% print(gcf, full_path, '-djpeg', '-r300');

%%
% figure

y_data = averageEdgeLengths;
x_data = scan_parameter_values;

% 分组计算
unique_x = unique(x_data);
means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);

% 绘制
% figure;
subplot(2, 2, 4);
errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'MarkerFaceColor', 'blue');
% xlabel('B (G)');
xlabel('\alpha (deg)');
ylabel('Average Spacing');
% title('分组统计图');
grid on;
theme("light")
set(gca, 'FontSize', 16);
xlim([min(x_data) max(x_data)])

%%
function [aspect_ratio] = calculate_aspect_ratios_simple(CC_split)
    
    % 获取区域属性，包括方向
    stats = regionprops(CC_split, 'BoundingBox', 'Centroid', 'Area', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
    
    num_regions = length(stats);
    
    aspect_ratio = zeros(num_regions, 1);
    
    for i = 1:num_regions
        major_length = stats(i).MajorAxisLength;
        minor_length = stats(i).MinorAxisLength;
        aspect_ratio(i) = major_length / minor_length;     
    end   
end


%%
% % 假设已有数据:
% % averageSumOD - 大小为(n,m)的二维数组
% % scan_parameter_values - 长度为n的一维数组
% 
% % 步骤1: 对行按照scan_parameter_values排序
% % 首先获取唯一的scan_parameter_values和对应的索引
% [unique_scan_values, ~, idx] = unique(scan_parameter_values, 'stable');
% 
% % 计算每个唯一值对应的平均索引位置
% group_indices = accumarray(idx, (1:length(scan_parameter_values))', [], @(x) {x});
% 
% % 创建新的排序索引
% sorted_row_indices = [];
% for i = 1:length(group_indices)
%     sorted_row_indices = [sorted_row_indices; group_indices{i}];
% end
% 
% % 按照排序后的索引重新排列数组
% sorted_averageSumOD = detectedArea(sorted_row_indices, :);
% sorted_scan_values = scan_parameter_values(sorted_row_indices);
% 
% % 步骤2: 对每列按照从大到小排序
% % 对每一列单独排序（降序）
% [~, col_indices] = sort(sorted_averageSumOD, 1, 'descend');
% 
% % 由于列排序会影响行数据的一致性，这里我们选择对转置后的矩阵按行排序
% % 这样能保持每行数据的完整性
% sorted_averageSumOD_final = zeros(size(sorted_averageSumOD));
% for col = 1:size(sorted_averageSumOD, 2)
%     sorted_averageSumOD_final(:, col) = sorted_averageSumOD(col_indices(:, col), col);
% end
% 
% % 步骤3: 绘制线图
% figure;
% hold on;
% 
% % 为不同的曲线设置不同的颜色
% colors = lines(size(sorted_averageSumOD_final, 1));
% 
% % 绘制每条曲线
% for i = 1:size(sorted_averageSumOD_final, 1)
%     plot(sorted_averageSumOD_final(i, :), ...
%          'Color', colors(i, :), ...
%          'LineWidth', 1.5, ...
%          'DisplayName', sprintf('B = %.2f', sorted_scan_values(i)));
% end
% 
% % 添加图形标签和标题
% xlabel('Detected Area Index');
% ylabel('Average OD in Detected Area');
% % title('按扫描参数排序的 averageSumOD 曲线');
% legend('show', 'Location', 'best');
% grid on;
% 
% % 美化图形
% set(gca, 'FontSize', 12);
% hold off;

%%
% function visualize_rotated_bounding_boxes(CC_split, original_image)
% % 在原图上可视化每个区域的旋转最小外接矩形
% % 这需要计算每个区域的方向
% 
%     % 准备可视化图像
%     if islogical(original_image)
%         vis_image = im2uint8(original_image);
%         vis_image = cat(3, vis_image, vis_image, vis_image);
%     else
%         if size(original_image, 3) == 1
%             vis_image = cat(3, original_image, original_image, original_image);
%         else
%             vis_image = original_image;
%         end
%     end
% 
%     % 获取区域属性，包括方向
%     stats = regionprops(CC_split, 'BoundingBox', 'Centroid', 'Area', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
% 
%     num_regions = length(stats);
% 
%     figure;
%     imshow(vis_image);
%     hold on;
% 
%     colors = lines(num_regions);
% 
%     for i = 1:num_regions
%         centroid = stats(i).Centroid;
%         orientation = stats(i).Orientation;
%         major_length = stats(i).MajorAxisLength;
%         minor_length = stats(i).MinorAxisLength;
%         area = stats(i).Area;
% 
%         % 计算长宽比
%         aspect_ratio = major_length / minor_length;
% 
%         % 计算旋转矩形的四个顶点
%         angle = deg2rad(orientation); % 转换为弧度，注意方向
%         cos_angle = cos(angle);
%         sin_angle = sin(angle);
% 
%         half_length = major_length / 2;
%         half_width = minor_length / 2;
% 
%         % 四个顶点的相对坐标
%         corners_relative = [
%             -half_length, -half_width;
%             half_length, -half_width;
%             half_length, half_width;
%             -half_length, half_width
%         ];
% 
%         % 旋转顶点
%         rotation_matrix = [cos_angle, -sin_angle; sin_angle, cos_angle];
%         corners_rotated = corners_relative * rotation_matrix;
% 
%         % 平移至中心位置
%         corners_absolute = corners_rotated + centroid;
% 
%         % 绘制旋转矩形
%         patch(corners_absolute(:,1), corners_absolute(:,2), 'r', ...
%             'FaceColor', 'none', ...
%             'EdgeColor', colors(i,:), ...
%             'LineWidth', 2);
% 
%         % 绘制中心点和主轴
%         plot(centroid(1), centroid(2), 'o', ...
%             'MarkerSize', 6, ...
%             'MarkerFaceColor', colors(i,:), ...
%             'MarkerEdgeColor', 'white');
% 
%         % 绘制主轴方向
%         end_point = centroid + [cos_angle * half_length, -sin_angle * half_length];
%         plot([centroid(1), end_point(1)], [centroid(2), end_point(2)], ...
%             'Color', colors(i,:), 'LineWidth', 2);
% 
%     end   
%     hold off;
% end
% 
% % 使用示例
% % 生成示例图像或读取你的图像
% BW = false(size(imgCropped));
% for i=1:CC_split.NumObjects
%     BW(cell2mat(CC_split.PixelIdxList(i))) = true;
% end
% 
% % 方法2：绘制旋转的最小外接矩形
% visualize_rotated_bounding_boxes(CC_split, BW);

%%
function [allRotatedPoints, edgeLengths, edgesPlot] = simpleAlignedNewtonDiagram(points)

    % 获取内层点
    DT = delaunayTriangulation(points);

    if isempty(DT.Points) || length(DT.Points) <3
        allRotatedPoints = [];
        edgeLengths = [];
        edgesPlot = [];
        return
    end

    convexHullIdx = convexHull(DT);
    isOuterPoint = false(size(points, 1), 1);
    isOuterPoint(convexHullIdx) = true;
    innerPoints = points(~isOuterPoint, :);
    innerPointIndices = find(~isOuterPoint);

    edges = DT.edges;
    edgesPlot = [];
    allAngles = [];  % 存储所有夹角
    
    if isempty(innerPoints)

        % 当innerPoints为空时，为每个点保留最短的边
        % 预分配存储空间
        numPoints = size(points, 1);
        shortestEdgeInfo = cell(numPoints, 1); % 存储每个点的最短边信息
        edgeLengths = [];
        edgeInner = [];
        
        % 第一遍：找到每个点的最短边
        for i = 1:size(edges, 1)
            edge = edges(i, :);
            point1 = edge(1);
            point2 = edge(2);
            p1 = points(point1, :);
            p2 = points(point2, :);
            currentLength = norm(p1 - p2);
            
            % 更新点1的最短边
            if isempty(shortestEdgeInfo{point1}) || ...
               currentLength < shortestEdgeInfo{point1}.length
                shortestEdgeInfo{point1} = struct(...
                    'edge', edge, ...
                    'length', currentLength);
            end
            
            % 更新点2的最短边
            if isempty(shortestEdgeInfo{point2}) || ...
               currentLength < shortestEdgeInfo{point2}.length
                shortestEdgeInfo{point2} = struct(...
                    'edge', edge, ...
                    'length', currentLength);
            end
        end
        
        % 第二遍：收集所有最短边（去重）
        collectedEdges = []; % 存储去重后的边索引
        for i = 1:numPoints
            if ~isempty(shortestEdgeInfo{i})
                currentEdge = sort(shortestEdgeInfo{i}.edge);
                
                % 检查是否已经收集过这条边
                isNew = true;
                for j = 1:size(collectedEdges, 1)
                    if all(collectedEdges(j, :) == currentEdge)
                        isNew = false;
                        break;
                    end
                end
                
                if isNew
                    collectedEdges = [collectedEdges; currentEdge];
                    edgeLengths = [edgeLengths; shortestEdgeInfo{i}.length];
                    edgeInner = [edgeInner; shortestEdgeInfo{i}.edge];
                end
            end
        end

    else

        edgeLengths = [];
        edgeInner = [];

        for i = 1:size(edges, 1)
            edge = edges(i, :);
            p1 = points(edge(1), :);
            p2 = points(edge(2), :);

            % 如果这条边不是连接两个外层点的，则保留
            if ~(isOuterPoint(edge(1)) && isOuterPoint(edge(2)))
                edgeLengths = [edgeLengths; norm(p1 - p2)];
                edgeInner = [edgeInner; edge(1), edge(2)];
            end
        end

    end

    numInner = size(innerPoints, 1);
    
    % 存储旋转后的点
    allRotatedPoints = [];
    
    % 处理每个内层点
    for i = 1:numInner
        centerPointIdx = innerPointIndices(i);
        centerPoint = points(centerPointIdx, :);
        
        % 找到与当前点相连的边
        connectedEdges = edgeInner(any(edgeInner == centerPointIdx, 2), :);
        
        % 获取所有邻居点索引（排除自身）
        neighborIndices = unique(connectedEdges(connectedEdges ~= centerPointIdx));
        
        if isempty(neighborIndices)
        %     fprintf('内层点 %d 没有相连的邻居\n', i);
            continue;
        end

        % 计算到每个邻居的距离
        neighborPoints = points(neighborIndices, :);
        distances = sqrt(sum((neighborPoints - centerPoint).^2, 2));

        % 找到最近的邻居（第一条边）
        [minDist, minIdx] = min(distances);
        nearestNeighborIdx = neighborIndices(minIdx);
        nearestPoint = neighborPoints(minIdx, :);
        neighborPoints(minIdx, :) = [];
        
        % 计算旋转角度（使最近邻点位于0度方向）
        vectorToNearest = nearestPoint - centerPoint;
        rotationAngle = -atan2(vectorToNearest(2), vectorToNearest(1));
        
        % 旋转所有邻居点
        rotationMatrix = [cos(rotationAngle), -sin(rotationAngle);
                         sin(rotationAngle), cos(rotationAngle)];
        
        % 计算邻居点相对于中心点的坐标
        neighborVectors = neighborPoints - centerPoint;
        rotatedNeighbors = (rotationMatrix * neighborVectors')';
        % rotatedNeighbors = rotatedNeighbors(abs(rotatedNeighbors(:, 2)) > 1e-3, :);
        % 存储旋转后的点
        allRotatedPoints = [allRotatedPoints; rotatedNeighbors];
    end

    edgesPlot = edgeInner;

    % % 创建图形
    % figure('Position', [100, 100, 800, 800]);
    % 
    % if ~isempty(allRotatedPoints)
    %     % 转换为极坐标
    %     [theta, rho] = cart2pol(allRotatedPoints(:,1), allRotatedPoints(:,2));
    % 
    %     % 创建极坐标图
    %     polarscatter(theta, rho, 30, 'filled', 'MarkerFaceAlpha', 0.6);
    %     hold on;
    % 
    %     % 添加平均半径圆
    %     meanRho = mean(rho);
    %     theta_circle = linspace(0, 2*pi, 100);
    %     polarplot(theta_circle, meanRho*ones(size(theta_circle)), ...
    %              'r-', 'LineWidth', 2);
    % 
    %     % 添加标准差范围
    %     stdRho = std(rho);
    %     polarplot(theta_circle, (meanRho+stdRho)*ones(size(theta_circle)), ...
    %              'g--', 'LineWidth', 1);
    %     polarplot(theta_circle, max(0, meanRho-stdRho)*ones(size(theta_circle)), ...
    %              'g--', 'LineWidth', 1);
    % 
    %     title(sprintf('点 %d\n半径: %.2f±%.2f', i, meanRho, stdRho));
    % else
    %     title(sprintf('点 %d\n无相邻点', i));
    % end
    % 
    % sgtitle('内层点的简化牛顿图');
end

% 使用示例： 
% simpleAlignedNewtonDiagram(centers);

%%
function [allAngles, edgeLengths, edgesPlot] = edgeAngleHistogram(points)
% edgeAngleHistogram - 计算并显示被选中边之间夹角的直方图
%
% 输入:
%   points - N×2矩阵，包含点的坐标
%
% 输出:
%   allAngles - 所有计算得到的夹角（度数）
%   edgeLengths - 被选中边的长度
%   edgesPlot - 被选中的边索引

    % 获取内层点
    DT = delaunayTriangulation(points);

    if isempty(DT.Points) || length(DT.Points) < 3
        allAngles = [];
        edgeLengths = [];
        edgesPlot = [];
        fprintf('点数不足，无法进行三角剖分\n');
        return
    end

    convexHullIdx = convexHull(DT);
    isOuterPoint = false(size(points, 1), 1);
    isOuterPoint(convexHullIdx) = true;
    innerPoints = points(~isOuterPoint, :);
    innerPointIndices = find(~isOuterPoint);

    edges = DT.edges;
    edgesPlot = [];
    allAngles = [];  % 存储所有夹角

    if isempty(innerPoints)
        % % 当innerPoints为空时，为每个点保留最短的边
        % numPoints = size(points, 1);
        % shortestEdgeInfo = cell(numPoints, 1);
        % edgeLengths = [];
        % edgeInner = [];
        % 
        % % 第一遍：找到每个点的最短边
        % for i = 1:size(edges, 1)
        %     edge = edges(i, :);
        %     point1 = edge(1);
        %     point2 = edge(2);
        %     p1 = points(point1, :);
        %     p2 = points(point2, :);
        %     currentLength = norm(p1 - p2);
        % 
        %     % 更新点1的最短边
        %     if isempty(shortestEdgeInfo{point1}) || ...
        %        currentLength < shortestEdgeInfo{point1}.length
        %         shortestEdgeInfo{point1} = struct(...
        %             'edge', edge, ...
        %             'length', currentLength);
        %     end
        % 
        %     % 更新点2的最短边
        %     if isempty(shortestEdgeInfo{point2}) || ...
        %        currentLength < shortestEdgeInfo{point2}.length
        %         shortestEdgeInfo{point2} = struct(...
        %             'edge', edge, ...
        %             'length', currentLength);
        %     end
        % end
        % 
        % % 第二遍：收集所有最短边（去重）
        % collectedEdges = [];
        % for i = 1:numPoints
        %     if ~isempty(shortestEdgeInfo{i})
        %         currentEdge = sort(shortestEdgeInfo{i}.edge);
        % 
        %         % 检查是否已经收集过这条边
        %         isNew = true;
        %         for j = 1:size(collectedEdges, 1)
        %             if all(collectedEdges(j, :) == currentEdge)
        %                 isNew = false;
        %                 break;
        %             end
        %         end
        % 
        %         if isNew
        %             collectedEdges = [collectedEdges; currentEdge];
        %             edgeLengths = [edgeLengths; shortestEdgeInfo{i}.length];
        %             edgeInner = [edgeInner; shortestEdgeInfo{i}.edge];
        %         end
        %     end
        % end

        edgeLengths = [];
        edgeInner = [];

        for i = 1:size(edges, 1)
            edge = edges(i, :);
            p1 = points(edge(1), :);
            p2 = points(edge(2), :);

            edgeLengths = [edgeLengths; norm(p1 - p2)];
            edgeInner = [edgeInner; edge(1), edge(2)];
        end

    else
        edgeLengths = [];
        edgeInner = [];

        for i = 1:size(edges, 1)
            edge = edges(i, :);
            p1 = points(edge(1), :);
            p2 = points(edge(2), :);

            % 如果这条边不是连接两个外层点的，则保留
            if ~(isOuterPoint(edge(1)) && isOuterPoint(edge(2)))
                edgeLengths = [edgeLengths; norm(p1 - p2)];
                edgeInner = [edgeInner; edge(1), edge(2)];
            end
        end
    end

    % 计算每条边的向量表示
    vectors = zeros(size(edgeInner, 1), 2);
    for i = 1:size(edgeInner, 1)
        p1 = points(edgeInner(i, 1), :);
        p2 = points(edgeInner(i, 2), :);
        vectors(i, :) = p2 - p1;  % 边向量
    end

    % 计算所有边对之间的夹角
    numEdges = size(vectors, 1);
    
    if numEdges < 2
        fprintf('边数不足，无法计算夹角\n');
        return;
    end

    % 计算夹角
    for i = 1:numEdges-1
        for j = i+1:numEdges
            % 获取两个向量
            v1 = vectors(i, :);
            v2 = vectors(j, :);
            
            % 归一化向量
            norm_v1 = norm(v1);
            norm_v2 = norm(v2);
            
            if norm_v1 > 0 && norm_v2 > 0
                % 计算点积
                dot_product = dot(v1/norm_v1, v2/norm_v2);
                
                % 限制点积范围，防止数值误差
                dot_product = max(min(dot_product, 1), -1);
                
                % 计算夹角（弧度）
                angle_rad = acos(dot_product);
                
                % 转换为度数
                angle_deg = angle_rad * 180 / pi;
                
                % 添加到夹角列表（可以选择保留锐角或所有角度）
                % 保留锐角（0-90度）
                % if angle_deg > 90
                %     angle_deg = 180 - angle_deg;  % 转换为锐角
                % end
                
                allAngles = [allAngles; angle_deg];
            end
        end
    end

    edgesPlot = edgeInner;

end

edgeAngleHistogram(centers)

%%
function plotSeparatedAngleDistributions(allRotatedPoints, scan_parameter_values, numBins)
    % 将每个扫描参数的角度分布绘制在不同的子图中
    % allRotatedPoints: 三维数组 (scan_params × points × 2)
    % scan_parameter_values: 扫描参数值（可能有重复）
    % numBins: 角度区间的数量（扇区数），默认36（每10度一个扇区）
    
    if nargin < 3
        numBins = 36; % 默认36个扇区，每10度一个
    end
    
    % 1. 获取唯一扫描参数值
    [unique_scans, ~, idx_to_unique] = unique(scan_parameter_values);
    n_unique_scans = length(unique_scans);
    
    fprintf('发现 %d 个唯一扫描参数（原始 %d 个参数）\n', ...
            n_unique_scans, length(scan_parameter_values));
    
    % 2. 计算子图的行列布局
    % 计算合适的子图布局
    if n_unique_scans <= 4
        nRows = 1;
        nCols = n_unique_scans;
    elseif n_unique_scans <= 9
        nRows = ceil(sqrt(n_unique_scans));
        nCols = ceil(n_unique_scans / nRows);
    else
        nRows = ceil(n_unique_scans / 5);
        nCols = 5;
    end
    
    % 3. 创建图形
    figure('Position', [100, 100, 1600, 900]);
    
    % 4. 生成颜色映射
    if n_unique_scans <= 10
        % 对于少数参数，使用预定义的鲜艳颜色
        colors = lines(n_unique_scans);
    else
        % 对于多个参数，使用jet或parula颜色映射
        colors = jet(n_unique_scans);
    end
    
    % 5. 为每个唯一参数创建子图
    for i = 1:n_unique_scans
        % 选择当前子图
        ax = subplot(nRows, nCols, i, polaraxes);
        
        % 找到所有等于当前唯一值的索引
        same_param_idx = find(idx_to_unique == i);
        
        % 合并这些索引对应的所有点
        temp_points = [];
        for j = 1:length(same_param_idx)
            idx = same_param_idx(j);
            current_points = squeeze(allRotatedPoints(idx, :, :));
            
            % 移除无效数据
            valid_idx = all(isfinite(current_points), 2);
            valid_idx(all(current_points == 0, 2)) = false;
            current_points = current_points(valid_idx, :);
            
            % 添加到临时存储
            temp_points = [temp_points; current_points];
        end
        
        % 如果有点数据，绘制分布
        if ~isempty(temp_points)
            % 分离x和y坐标
            x_data = temp_points(:, 1);
            y_data = temp_points(:, 2);
            
            % 计算每个点的角度
            angles_rad = atan2(y_data, x_data);
            angles_deg = mod(rad2deg(angles_rad), 360);
            
            % 计算角度直方图（计数）
            angle_edges = linspace(0, 360, numBins + 1);
            angle_centers = angle_edges(1:end-1) + (angle_edges(2) - angle_edges(1))/2;
            counts = histcounts(angles_deg, angle_edges);
            
            % % 计算百分比
            total_points = length(angles_deg);
            % percentages = counts / total_points * 100;
            
            % 绘制极坐标图
            % 将角度转换为弧度（用于极坐标）
            theta_rad = deg2rad(angle_centers);
            
            % 闭合曲线（连接首尾）
            % theta_plot = [theta_rad, 2*pi];
            % r_plot = [counts, counts(1)];

            counts = counts + fliplr(counts);
            theta_plot = theta_rad;
            r_plot = counts;
            
            % 创建极坐标图
            polarplot(ax, theta_plot, r_plot, '-', ...
                     'Color', colors(i, :), ...
                     'LineWidth', 2);
            hold on;
            
            % 设置极坐标图属性
            ax.ThetaZeroLocation = 'right';  % 0度在右侧
            ax.ThetaDir = 'counterclockwise'; % 逆时针方向
            ax.RGrid = 'on';
            ax.GridLineStyle = '--';
            ax.GridAlpha = 0.3;
            ax.ThetaTick = 0:30:330;
            ax.ThetaTickLabel = {'0°', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°', '270°', '300°', '330°'};
            
            % 添加标题和信息
            title_str = sprintf('B = %.4f(G), n = %d', unique_scans(i), total_points);
            title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
            
            fprintf('参数 %d (值=%.4f): %d 个点\n', ...
                    i, unique_scans(i), total_points);
        else
            % 创建空极坐标图
            title(sprintf('B = %.4f (G)\nNo Data', unique_scans(i)), ...
                  'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
            text(0.5, 0.5, 'No Data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
        end

        thetalim([0 180]);  % 只显示0°到180°
        
        % 设置当前子图边框颜色
        set(gca, 'LineWidth', 1.5, 'Box', 'on');
    end

    theme("light")
    
    fprintf('\n绘图完成: %d 个子图已创建\n', n_unique_scans);
end

plotSeparatedAngleDistributions(allRotatedPoints, scan_parameter_values, 36);

%%
function plotAngleHistograms(allAngles, scan_parameter_values, numBins)
    % 将每个扫描参数的角度分布绘制在不同的子图中（直角坐标系）
    % allAngles: 二维数组 (scan_params × data_points)，包含角度值（度）
    % scan_parameter_values: 扫描参数值（可能有重复）
    % numBins: 直方图的箱数，默认18（每10度一个箱）
    
    if nargin < 3
        numBins = 18; % 默认18个箱，每10度一个
    end
    
    % 1. 验证输入数据维度
    n_scans = length(scan_parameter_values);
    if size(allAngles, 1) ~= n_scans
        error('allAngles的第一个维度大小(%d)与scan_parameter_values长度(%d)不匹配', ...
              size(allAngles, 1), n_scans);
    end
    
    % 获取数组维度信息
    [n_scans, n_points] = size(allAngles);
    fprintf('数组维度: %d × %d\n', n_scans, n_points);
    
    % 2. 获取唯一扫描参数值
    [unique_scans, ~, idx_to_unique] = unique(scan_parameter_values);
    n_unique_scans = length(unique_scans);
    
    fprintf('发现 %d 个唯一扫描参数（原始 %d 个参数）\n', ...
            n_unique_scans, length(scan_parameter_values));
    
    % 3. 计算子图的行列布局
    % 计算合适的子图布局
    if n_unique_scans <= 4
        nRows = 1;
        nCols = n_unique_scans;
    elseif n_unique_scans <= 9
        nRows = ceil(sqrt(n_unique_scans));
        nCols = ceil(n_unique_scans / nRows);
    else
        nRows = ceil(n_unique_scans / 5);
        nCols = 5;
    end
    
    % 4. 创建图形
    figure('Position', [100, 100, 1600, 900]);
    
    % 5. 生成颜色映射
    if n_unique_scans <= 10
        % 对于少数参数，使用预定义的鲜艳颜色
        colors = lines(n_unique_scans);
    else
        % 对于多个参数，使用jet或parula颜色映射
        colors = jet(n_unique_scans);
    end
    
    % 6. 为每个唯一参数创建子图
    for i = 1:n_unique_scans
        % 选择当前子图
        ax = subplot(nRows, nCols, i);
        
        % 找到所有等于当前唯一值的索引
        same_param_idx = find(idx_to_unique == i);
        
        % 合并这些索引对应的所有角度数据（从二维数组中提取）
        temp_angles = [];
        for j = 1:length(same_param_idx)
            idx = same_param_idx(j);
            
            % 提取当前扫描参数的所有角度数据（第二维）
            current_angles = allAngles(idx, :);
            
            % 移除无效数据 (NaN, Inf, 超出合理范围的值，以及0值)
            % 注意：这里排除了0值，因为您要求不统计0值
            valid_idx = isfinite(current_angles) & ...
                       (current_angles > 0) & ...     % 排除0值
                       (current_angles <= 360);       % 只保留0-360度的值
            valid_angles = current_angles(valid_idx);
            
            % 添加到总的角度数据
            temp_angles = [temp_angles, valid_angles];
        end
        
        % 如果有角度数据，绘制直方图
        if ~isempty(temp_angles)
            % 将角度归一化到0-360度范围
            angles_deg = mod(temp_angles, 360);
            
            % 确保归一化后没有0值（如果原始值接近360，mod后可能为0）
            angles_deg(angles_deg == 0) = 360;
            
            % 映射到0-180度范围（对称处理）
            mapped_angles = angles_deg;
            mapped_angles(mapped_angles > 180) = 360 - mapped_angles(mapped_angles > 180);
            
            % 计算直方图（0-180度范围）
            [counts, bin_edges] = histcounts(mapped_angles, numBins, 'BinEdges', 0:180/numBins:180);
            bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

            counts = counts/max(counts);
            
            % 绘制条形图
            bar(ax, bin_centers, counts, 'FaceColor', colors(i, :), 'EdgeColor', colors(i, :), 'FaceAlpha', 0.7);
            hold on;
            
            % 添加平滑曲线（可选）
            smooth_x = linspace(0, 180, 500);
            smooth_counts = interp1([bin_centers(1)-180/numBins, bin_centers, bin_centers(end)+180/numBins], ...
                                   [counts(end), counts, counts(1)], smooth_x, 'pchip');
            plot(ax, smooth_x, smooth_counts, 'Color', colors(i, :), 'LineWidth', 2);
            
            % 设置坐标轴属性
            xlim([0, 180]);
            xlabel('Degree', 'FontSize', 10);
            ylabel('Count', 'FontSize', 10);
            
            % 设置x轴刻度
            xticks(0:30:180);
            xticklabels({'0°', '30°', '60°', '90°', '120°', '150°', '180°'});
            
            % 添加网格
            grid on;
            grid minor;
            
            % 添加标题和信息
            title_str = sprintf('B = %.4f (G), n = %d', unique_scans(i), length(mapped_angles));
            title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
            
            fprintf('参数 %d (值=%.4f): %d 个角度值 (排除0值后)\n', ...
                    i, unique_scans(i), length(mapped_angles));
            
            % 计算并显示统计信息
            mean_angle = mean(mapped_angles);
            std_angle = std(mapped_angles);
            median_angle = median(mapped_angles);
            
            % 在图中添加统计信息
            text_x = 0.98;
            text_y = 0.95;
            % annotation_text = sprintf('均值: %.1f°\n中位数: %.1f°\n标准差: %.1f°', ...
            %                           mean_angle, median_angle, std_angle);
            % text(ax, text_x*180, text_y*max(counts), annotation_text, ...
            %      'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
            %      'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 2);

            % ylim([0, 70]);
        else
            % 创建空图
            title(sprintf('α = %.4f\nNo Data', unique_scans(i)), ...
                  'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
            text(0.5, 0.5, 'No Data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
            xlim([0, 180]);
            ylim([0, 1]);
            xlabel('Degree', 'FontSize', 10);
            ylabel('count', 'FontSize', 10);
            xticks(0:30:180);
        end
        
        % 设置当前子图边框颜色
        set(gca, 'LineWidth', 1, 'Box', 'on');
    end

    % 为整个图形添加标题
    sgtitle('Angular distribution histogram', 'FontSize', 14, 'FontWeight', 'bold');
    
    theme("light")
    
    fprintf('\n绘图完成: %d 个子图已创建\n', n_unique_scans);
end

% 调用函数
plotAngleHistograms(allAngles, scan_parameter_values, 20);