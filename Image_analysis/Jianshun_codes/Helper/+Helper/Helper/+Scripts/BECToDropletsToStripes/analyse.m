folderPath = "C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2025\10\31\0002\AnalysisPlot\Data\";
file_list = dir(fullfile(folderPath, '*.mat'));

load(fullfile(folderPath, file_list(1).name));

number = zeros(length(scan_parameter_values), 1);
area = zeros(length(scan_parameter_values), 1);
sumOD = zeros(length(scan_parameter_values), 1);

for i = 1:16%length(file_list)
    load(fullfile(folderPath, file_list(i).name));
    number(i) = CC_split.NumObjects;

    if CC_split.NumObjects == 0
        area(i) = 0;
    else
        % 计算每个区域的面积（像素数量）
        areas = cellfun(@numel, CC_split.PixelIdxList);
    
        % 计算平均面积
        area(i) = mean(areas);
    
    end

    numRegions = length(CC_split.PixelIdxList);

    grayValueSums = zeros(numRegions, 1);

    % 遍历每个区域，计算灰度值总和
    for k = 1:numRegions
        % 获取当前区域的像素索引
        pixelIdx = cell2mat(CC_split.PixelIdxList(k));
    
        % 计算该区域的灰度值总和
        grayValueSums(k) = sum(imgCropped(pixelIdx));
    end

    sumOD(i) = mean(grayValueSums);

    % img = double(OD);
    % calculateWeightedCenters(CC_split, imgCropped)
end

%%
y_data = sumOD;
x_data = scan_parameter_values;

% 分组计算
unique_x = unique(x_data);
means = arrayfun(@(x) mean(y_data(x_data == x)), unique_x);
stds = arrayfun(@(x) std(y_data(x_data == x)), unique_x);

% 绘制
figure;
errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'MarkerFaceColor', 'blue');
xlabel('B (G)');
ylabel('Detected Area Size');
% title('分组统计图');
grid on;

%%
y_data = number;
x_data = scan_parameter_values;

% 分组计算
unique_x = unique(x_data);
means = arrayfun(@(x) mean(y_data(x_data == x)), unique_x);
stds = arrayfun(@(x) std(y_data(x_data == x)), unique_x);

% 绘制
figure;
errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'MarkerFaceColor', 'blue');
xlabel('B (G)');
ylabel('Detected Area Number');
% title('分组统计图');
grid on;

%%
figure

plot(scan_parameter_values, number, 'o')

%%

centers = calculateWeightedCenters(CC_split, imgCropped);

DT = delaunayTriangulation(centers);

%%

% 创建标签图像
labeledImage = labelmatrix(CC_split);

% 可视化：在一张图上标出中心
figure;
imshow(imgCropped, []);
hold on;

% 显示彩色标签图像
imshow(label2rgb(labeledImage, 'jet', 'w', 'shuffle'));
hold on;

% 绘制区域边界
boundaries = bwboundaries(imgCropped);
for k = 1:length(boundaries)
    boundary = boundaries{k};
    plot(boundary(:, 2), boundary(:, 1), 'k', 'LineWidth', 1);
end

% 绘制Delaunay三角网格
triplot(DT, 'g-', 'LineWidth', 1.5);

% 标记中心点
plot(centers(:, 1), centers(:, 2), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
plot(centers(:, 1), centers(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% 添加区域编号标签
for i = 1:size(centers, 1)
    text(centers(i, 1) + 5, centers(i, 2) - 5, num2str(i), ...
         'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold');
end

title('连通区域加权中心', 'FontSize', 14);
hold off;

% 可选：在命令窗口显示坐标
fprintf('各区域加权中心坐标（列，行）：\n');
for i = 1:size(centers, 1)
    fprintf('区域 %d: (%.2f, %.2f)\n', i, centers(i, 1), centers(i, 2));
end

%%

folderPath = "C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2025\10\31\0000\AnalysisPlot\Data\";
file_list = dir(fullfile(folderPath, '*.mat'));

number = zeros(length(scan_reference_values));

for i = 1:1%length(file_list)
    load(fullfile(folderPath, file_list(i).name));
    number(i) = CC_split.NumObjects;
end

figure

plot(scan_reference_values, number')

%%

figure;
imshow(imgCropped)
colormap jet