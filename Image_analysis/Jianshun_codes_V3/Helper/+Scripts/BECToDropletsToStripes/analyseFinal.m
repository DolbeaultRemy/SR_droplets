% main_analysis.m
% 主要分析脚本：加载数据、处理图像区域、计算统计数据并绘图

clear all; close all; clc;

%% 1. 数据加载和初始化
% 设置数据文件夹路径

global folderPath file_list scan_parameter_values number ...
           averageArea averageSumOD sumOD detectedArea averageOD ...
           selectFlag aspectRatio averageAspectRatio allRotatedPoints ...
           edgeLengths averageEdgeLengths allAngles plotFlag...
           CC_split imgCropped dataset_label dataset_description...
           ResultPath

plotFlag = false;  % 控制是否绘制每个文件的详细图像

dataFolder = [];

% ResultPath = "D:\Jianshun\Data\2026\01\01\0001\AnalysisPlot\";
% dataset_label = "180";
% dataset_description = "202601010001, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\01\0002\AnalysisPlot\";
% dataset_label = "180";
% dataset_description = "202601010002, MS to BEC";
% main()
% dataFolder = [dataFolder, ResultPath];

ResultPath = "C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2026\01\01\0005\AnalysisPlot\";
dataset_label = "165";
dataset_description = "202601010005, BEC to MS";
main()
dataFolder = [dataFolder, ResultPath];

% ResultPath = "D:\Jianshun\Data\2026\01\02\0000\AnalysisPlot\";
% dataset_label = "165";
% dataset_description = "202601020000, MS to BEC";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\02\0004\AnalysisPlot\";
% dataset_label = "150";
% dataset_description = "202601020004, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\02\0012\AnalysisPlot\";
% dataset_label = "135";
% dataset_description = "202601020012, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\03\0003\AnalysisPlot\";
% dataset_label = "120";
% dataset_description = "202601030003, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\05\0004\AnalysisPlot\";
% dataset_label = "120";
% dataset_description = "202601050004, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\05\0005\AnalysisPlot\";
% dataset_label = "120";
% dataset_description = "202601050005, MS to BEC";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\06\0002\AnalysisPlot\";
% dataset_label = "110";
% dataset_description = "202601060002, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\11\0007\AnalysisPlot\";
% dataset_label = "140";
% dataset_description = "202601110007, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\12\0002\AnalysisPlot\";
% dataset_label = "130";
% dataset_description = "202601120002, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\12\0006\AnalysisPlot\";
% dataset_label = "125";
% dataset_description = "202601120006, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% ResultPath = "D:\Jianshun\Data\2026\01\13\0002\AnalysisPlot\";
% dataset_label = "140";
% dataset_description = "202601130002, BEC to MS";
% main()
% dataFolder = [dataFolder, ResultPath];
% 
% FinialResultPath = "D:\Jianshun\Data\2026\01\";
% combine_plot_data(dataFolder, FinialResultPath);

%% Main function
function main()

    global folderPath file_list scan_parameter_values number ...
           averageArea averageSumOD sumOD detectedArea averageOD ...
           selectFlag aspectRatio averageAspectRatio allRotatedPoints ...
           edgeLengths averageEdgeLengths allAngles plotFlag...
           CC_split imgCropped dataset_label dataset_description...
           ResultPath

    folderPath = ResultPath + "Data\";
    file_list = dir(fullfile(folderPath, '*.mat'));
    
    % 加载第一个文件以获取数据结构
    load(fullfile(folderPath, file_list(1).name));
    
    % 初始化参数
    plotIdx = length(file_list);
    
    % 预分配数据存储数组
    initialize_data_arrays();
    
    % 2. 主处理循环
    for i = 1:plotIdx
        process_single_file(i);
    end
    
    % 后处理：清除无效数据
    aspectRatio(~selectFlag) = 0;
    
    % 3. 绘制汇总统计图
    create_summary_plots();
end

%% 4. 辅助函数定义
function initialize_data_arrays()
    % 初始化所有数据存储数组
    
    % 全局变量声明（在函数内部使用persistent或作为主脚本变量）
    global number averageArea averageSumOD sumOD detectedArea ...
           averageOD selectFlag aspectRatio averageAspectRatio ...
           allRotatedPoints edgeLengths averageEdgeLengths allAngles ...
           scan_parameter_values
    
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
end

function process_single_file(file_idx)
    % 处理单个MAT文件
    
    global folderPath file_list scan_parameter_values number ...
           averageArea averageSumOD sumOD detectedArea averageOD ...
           selectFlag aspectRatio averageAspectRatio allRotatedPoints ...
           edgeLengths averageEdgeLengths allAngles plotFlag...
           CC_split imgCropped
    
    % 加载当前文件
    current_file = fullfile(folderPath, file_list(file_idx).name);
    load(current_file);
    
    % 处理分割区域
    if CC_split.NumObjects == 0
        averageArea(file_idx) = 0;
    else
        process_regions(file_idx);
    end
    
    % 筛选感兴趣的区域
    apply_selection_criteria(file_idx);
    
    % 计算区域统计数据
    calculate_region_statistics(file_idx);
    
    % 计算加权中心并分析几何特性
    if any(selectFlag(file_idx, :))
        centers = calculateWeightedCenters(CC_split, imgCropped, selectFlag(file_idx, :));
        analyze_geometric_properties(file_idx, centers);
    else
        centers = [];
    end
    
    % 可选：绘制当前文件的结果
    if plotFlag
        plot_single_file_results(file_idx, centers);
    end
end

function process_regions(file_idx)
    % 处理分割区域，计算面积和光密度
    
    global CC_split imgCropped detectedArea sumOD averageOD
    
    numRegions = length(CC_split.PixelIdxList);
    
    % 计算每个区域的面积（像素数量）
    areas = cellfun(@numel, CC_split.PixelIdxList);
    detectedArea(file_idx, 1:length(areas)) = areas;
    
    % 计算每个区域的光密度总和
    for k = 1:numRegions
        pixelIdx = cell2mat(CC_split.PixelIdxList(k));
        sumOD(file_idx, k) = sum(imgCropped(pixelIdx));
        averageOD(file_idx, k) = sumOD(file_idx, k) / detectedArea(file_idx, k);
    end
end

function apply_selection_criteria(file_idx)
    % 应用选择标准筛选感兴趣的区域
    
    global CC_split averageOD detectedArea selectFlag aspectRatio
    
    % 计算长宽比
    aspectRatioTemp = calculate_aspect_ratios_simple(CC_split);
    aspectRatio(file_idx, 1:length(aspectRatioTemp)) = aspectRatioTemp;
    
    % 基于光密度和面积筛选区域
    for j = 1:CC_split.NumObjects
        if averageOD(file_idx, j) > 2.5 && detectedArea(file_idx, j) > 30
            selectFlag(file_idx, j) = 1;
        end
    end
end

function calculate_region_statistics(file_idx)
    % 计算选定区域的统计数据
    
    global CC_split detectedArea sumOD aspectRatio selectFlag ...
           number averageArea averageSumOD averageAspectRatio
    
    % 初始化累加器
    total_area = 0;
    total_sumOD = 0;
    total_aspectRatio = 0;
    
    % 累加选定区域的数据
    for j = 1:CC_split.NumObjects
        if selectFlag(file_idx, j)
            total_area = total_area + detectedArea(file_idx, j);
            total_sumOD = total_sumOD + sumOD(file_idx, j);
            total_aspectRatio = total_aspectRatio + aspectRatio(file_idx, j);
        end
    end
    
    % 计算平均统计数据
    number(file_idx) = sum(selectFlag(file_idx, :));
    
    if number(file_idx) > 0
        averageArea(file_idx) = total_area / number(file_idx);
        averageSumOD(file_idx) = total_sumOD / number(file_idx);
        averageAspectRatio(file_idx) = total_aspectRatio / number(file_idx);
    else
        averageArea(file_idx) = 0;
        averageSumOD(file_idx) = 0;
        averageAspectRatio(file_idx) = 0;
    end
end

function analyze_geometric_properties(file_idx, centers)
    % 分析几何特性：牛顿图和角度直方图
    
    global allRotatedPoints edgeLengths averageEdgeLengths allAngles
    
    % 计算对齐的牛顿图
    [allRotatedPointSingle, edgeLength, ~] = simpleAlignedNewtonDiagram(centers);
    
    if ~isempty(allRotatedPointSingle)
        allRotatedPoints(file_idx, 1:size(allRotatedPointSingle, 1), :) = allRotatedPointSingle;
        edgeLengths(file_idx, 1:length(edgeLength)) = edgeLength;
        averageEdgeLengths(file_idx) = mean(edgeLength);
    end
    
    % 计算边角度直方图
    [allAngle, ~, ~] = edgeAngleHistogram(centers);
    if ~isempty(allAngle)
        allAngles(file_idx, 1:length(allAngle)) = allAngle;
    end
end

function plot_single_file_results(file_idx, centers)
    % 绘制单个文件的详细结果
    
    global CC_split imgCropped selectFlag scan_parameter_values folderPath
    
    % 创建选定区域的掩码
    mask = false(size(imgCropped));
    for k = 1:CC_split.NumObjects
        if selectFlag(file_idx, k)
            mask(CC_split.PixelIdxList{k}) = true;
        end
    end
    
    % 获取区域边界
    boundaries = bwboundaries(mask);
    
    % 创建图形
    figure(2);
    imagesc(imgCropped, [0, 3]);
    colormap jet;
    clim([-0.5, 3.5]);
    hold on;
    
    % 绘制区域边界
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2);
    end
    
    % 尝试计算Delaunay三角剖分
    try
        DT = delaunayTriangulation(centers);
        convexHullIndices = convexHull(DT);
        
        % 移除重复点
        if convexHullIndices(1) == convexHullIndices(end)
            convexHullIndices = convexHullIndices(1:end-1);
        end
        
        % 分离内部点和外部点
        isOuterPoint = false(size(centers, 1), 1);
        isOuterPoint(convexHullIndices) = true;
        innerPoints = centers(~isOuterPoint, :);
        outerPoints = centers(isOuterPoint, :);
        
        % 绘制Delaunay三角网格（蓝色）
        triplot(DT, 'b-', 'LineWidth', 2.5);
        
        % 计算选中的边（红色）
        plot_selected_edges(centers, DT, isOuterPoint);
        
        % 标记中心点
        plot(innerPoints(:, 1), innerPoints(:, 2), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
        plot(innerPoints(:, 1), innerPoints(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        plot(outerPoints(:, 1), outerPoints(:, 2), 'k+', 'MarkerSize', 15, 'LineWidth', 2);
        plot(outerPoints(:, 1), outerPoints(:, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    catch
        % 如果三角剖分失败，只绘制标题
    end
    
    % 设置图形属性
    title("B = " + sprintf('%.2f (G)', scan_parameter_values(file_idx)));
    theme("light");
    set(gca, 'FontSize', 16);
    set(gcf, 'Position', [-1200, 200, 900, 600]);
    hold off;
    
    % 保存图像
    save_plot_image(file_idx);
end

function plot_selected_edges(centers, DT, isOuterPoint)
    % 绘制选中的边为红色
    
    % 获取所有边
    edges = DT.edges;
    
    % 初始化选中的边
    selected_edges = [];
    
    % 获取内层点索引
    innerPointIndices = find(~isOuterPoint);
    
    % 如果存在内层点
    if ~isempty(innerPointIndices)
        % 计算每条边的长度
        edge_lengths = zeros(size(edges, 1), 1);
        for i = 1:size(edges, 1)
            edge = edges(i, :);
            p1 = centers(edge(1), :);
            p2 = centers(edge(2), :);
            edge_lengths(i) = norm(p1 - p2);
        end
        
        % 对每条边判断是否选中
        for i = 1:size(edges, 1)
            edge = edges(i, :);
            
            % 如果这条边不是连接两个外层点的，则选中
            if ~(isOuterPoint(edge(1)) && isOuterPoint(edge(2)))
                selected_edges = [selected_edges; edge];
            end
        end
    else
        % 如果没有内层点，为每个点保留最短的边
        numPoints = size(centers, 1);
        shortestEdgeInfo = cell(numPoints, 1);
        
        % 第一遍：找到每个点的最短边
        for i = 1:size(edges, 1)
            edge = edges(i, :);
            point1 = edge(1);
            point2 = edge(2);
            p1 = centers(point1, :);
            p2 = centers(point2, :);
            currentLength = norm(p1 - p2);
            
            % 更新点1的最短边
            if isempty(shortestEdgeInfo{point1}) || currentLength < shortestEdgeInfo{point1}.length
                shortestEdgeInfo{point1} = struct('edge', edge, 'length', currentLength);
            end
            
            % 更新点2的最短边
            if isempty(shortestEdgeInfo{point2}) || currentLength < shortestEdgeInfo{point2}.length
                shortestEdgeInfo{point2} = struct('edge', edge, 'length', currentLength);
            end
        end
        
        % 第二遍：收集所有最短边（去重）
        collectedEdges = [];
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
                    selected_edges = [selected_edges; shortestEdgeInfo{i}.edge];
                end
            end
        end
    end
    
    % 绘制选中的边为红色
    for i = 1:size(selected_edges, 1)
        edge = selected_edges(i, :);
        x_red = [centers(edge(1), 1), centers(edge(2), 1)];
        y_red = [centers(edge(1), 2), centers(edge(2), 2)];
        plot(x_red, y_red, 'r-', 'LineWidth', 3);
    end
end

function save_plot_image(file_idx)
    % 保存绘图图像到文件
    
    global folderPath scan_parameter_values
    
    save_path = fullfile(folderPath, "FinialPlot");
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    filename = sprintf('%d_B = %.2f (G).jpg', file_idx, scan_parameter_values(file_idx));
    full_path = fullfile(save_path, filename);
    print(gcf, full_path, '-djpeg', '-r300');
end

function create_summary_plots()
    % 创建汇总统计图
    
    global scan_parameter_values number averageArea averageAspectRatio ...
           averageEdgeLengths aspectRatio allRotatedPoints allAngles
    
    % 创建包含4个子图的图形窗口
    % figure('Position', [100, 100, 1200, 800]);
    figure
    
    % 子图1: 检测到的区域数量
    subplot(2, 2, 1);
    [number_means, number_stds, number_unique_x] = plot_grouped_statistics(...
        scan_parameter_values, number, 'α (deg)', 'Detected Area Number');
    
    % 子图2: 检测到的区域大小
    subplot(2, 2, 2);
    [area_means, area_stds, area_unique_x] = plot_grouped_statistics(...
        scan_parameter_values, averageArea, 'α (deg)', 'Detected Area Size');
    
    % 子图3: 长宽比分布
    subplot(2, 2, 3);
    [aspect_means, aspect_stds, aspect_unique_x] = plot_aspect_ratio_distribution();
    
    % 子图4: 平均间距
    subplot(2, 2, 4);
    [spacing_means, spacing_stds, spacing_unique_x] = plot_grouped_statistics(...
        scan_parameter_values, averageEdgeLengths, 'α (deg)', 'Average Spacing');
    
    % 保存汇总图
    save_result_plot("Sumarry");

    % 保存绘图数据到MAT文件
    save_plotting_data(number_means, number_stds, number_unique_x, ...
                       area_means, area_stds, area_unique_x, ...
                       aspect_means, aspect_stds, aspect_unique_x, ...
                       spacing_means, spacing_stds, spacing_unique_x);
    
    % 创建角度分布图
    plot_angle_distributions();
end

function save_plotting_data(number_means, number_stds, number_unique_x, ...
                            area_means, area_stds, area_unique_x, ...
                            aspect_means, aspect_stds, aspect_unique_x, ...
                            spacing_means, spacing_stds, spacing_unique_x)
    % 保存绘图数据到MAT文件
    
    global ResultPath dataset_label dataset_description
    
    % 创建数据结构
    plot_data = struct();

    % 添加数据标签
    % 可以根据需要修改这里的数据标签
    plot_data.dataset_label = dataset_label;  % 这里手动设置数据标签
    plot_data.description = dataset_description;  % 可选：添加描述
    
    % 区域数量数据
    plot_data.detected_number = struct();
    plot_data.detected_number.x = number_unique_x;
    plot_data.detected_number.mean = number_means;
    plot_data.detected_number.std = number_stds;
    plot_data.detected_number.xlabel = 'α (deg)';
    plot_data.detected_number.ylabel = 'Detected Area Number';
    
    % 区域面积数据
    plot_data.detected_area = struct();
    plot_data.detected_area.x = area_unique_x;
    plot_data.detected_area.mean = area_means;
    plot_data.detected_area.std = area_stds;
    plot_data.detected_area.xlabel = 'α (deg)';
    plot_data.detected_area.ylabel = 'Detected Area Size';
    
    % 长宽比数据
    plot_data.aspect_ratio = struct();
    plot_data.aspect_ratio.x = aspect_unique_x;
    plot_data.aspect_ratio.mean = aspect_means;
    plot_data.aspect_ratio.std = aspect_stds;
    plot_data.aspect_ratio.xlabel = 'α (deg)';
    plot_data.aspect_ratio.ylabel = 'Average Aspect Ratio';
    
    % 平均间距数据
    plot_data.average_spacing = struct();
    plot_data.average_spacing.x = spacing_unique_x;
    plot_data.average_spacing.mean = spacing_means;
    plot_data.average_spacing.std = spacing_stds;
    plot_data.average_spacing.xlabel = 'α (deg)';
    plot_data.average_spacing.ylabel = 'Average Spacing';
    
    % 保存扫描参数信息
    global scan_parameter_values
    plot_data.scan_parameters = struct();
    plot_data.scan_parameters.values = scan_parameter_values;
    plot_data.scan_parameters.xlabel = 'α (deg)';
    
    % 保存到文件
    save_filename = fullfile(ResultPath, 'plotting_data.mat');
    save(save_filename, 'plot_data');
    
    fprintf('绘图数据已保存到: %s\n', save_filename);
    fprintf('数据包含以下结构:\n');
    fprintf('  plot_data.detected_number - 检测到的区域数量\n');
    fprintf('  plot_data.detected_area - 检测到的区域大小\n');
    fprintf('  plot_data.aspect_ratio - 长宽比\n');
    fprintf('  plot_data.average_spacing - 平均间距\n');
    fprintf('  plot_data.scan_parameters - 扫描参数\n');
end

function [means, stds, unique_x] = plot_grouped_statistics(x_data, y_data, xlabel_str, ylabel_str)
    % 绘制分组统计数据（带误差棒）并返回统计值
    % 输入:
    %   x_data - x轴数据
    %   y_data - y轴数据
    %   xlabel_str - x轴标签
    %   ylabel_str - y轴标签
    % 输出:
    %   means - 每个x值的平均值
    %   stds - 每个x值的标准差
    %   unique_x - 唯一的x值
    
    % 分组计算
    unique_x = unique(x_data);
    means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
    stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);
    
    % 绘制误差棒图
    errorbar(unique_x, means, stds, 'o-', 'LineWidth', 2, ...
             'MarkerSize', 8, 'MarkerFaceColor', 'blue');
    
    % 设置坐标轴标签和属性
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    grid on;
    theme("light");
    set(gca, 'FontSize', 16);
    xlim([min(x_data), max(x_data)]);
end

function [means, stds, unique_x] = plot_aspect_ratio_distribution()
    % 绘制长宽比分布图（线图和伪色彩图使用同一坐标轴）并返回统计值
    % 输出:
    %   means - 每个x值的平均长宽比
    %   stds - 每个x值的标准差
    %   unique_x - 唯一的x值
    
    global scan_parameter_values aspectRatio averageAspectRatio
    
    % 1. 对aspectRatio的行以scan_parameter_values为关键字排序
    [sorted_scan, sort_idx] = sort(scan_parameter_values);
    sorted_aspectRatio = aspectRatio(sort_idx, :);
    
    % 2. 准备数据
    n = size(sorted_aspectRatio, 1);
    m = size(sorted_aspectRatio, 2);
    
    % 创建对应的x和y数据
    x_data = repelem(sorted_scan, m);
    y_data = sorted_aspectRatio';
    y_data = y_data(:);
    
    % 移除无效数据点
    valid_idx = isfinite(y_data) & isfinite(x_data);
    x_data = x_data(valid_idx);
    y_data = y_data(valid_idx);
    
    % 3. 设置y轴的bin数量
    y_bins = 50;
    
    % 4. 计算二维直方图
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
        y_values = y_data(x_data == current_x);
        
        if ~isempty(y_values)
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
    hold on;
    h = pcolor(X, Y, counts');
    set(h, 'EdgeColor', 'none');
    
    % 8. 设置colormap
    colormap(jet);
    colorbar;
    
    % 9. 在同一坐标轴上绘制平均长宽比曲线
    [means, stds] = plot_average_aspect_ratio_on_same_axis();
    
    % 10. 添加标签
    xlabel('α (deg)');
    ylabel('Aspect Ratio');
    title('Aspect Ratio Distribution vs Scan Parameter');
    
    hold off;
end

function [means, stds] = plot_average_aspect_ratio_on_same_axis()
    % 在同一坐标轴上绘制平均长宽比曲线并返回统计值
    % 输出:
    %   means - 每个x值的平均长宽比
    %   stds - 每个x值的标准差
    
    global scan_parameter_values averageAspectRatio
    
    y_data = averageAspectRatio;
    x_data = scan_parameter_values;
    
    % 分组计算
    unique_x = unique(x_data);
    means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
    stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);
    
    % 绘制误差棒图
    a = errorbar(unique_x, means, stds, 'ko-', 'LineWidth', 3, ...
                 'MarkerSize', 4, 'MarkerFaceColor', 'white', ...
                 "DisplayName", "Average value");
    
    % 添加图例
    legend(a, 'Location', 'best');
end

function plot_average_aspect_ratio_on_top()
    % 在长宽比分布图顶部绘制平均长宽比曲线
    
    global scan_parameter_values averageAspectRatio
    
    y_data = averageAspectRatio;
    x_data = scan_parameter_values;
    
    % 分组计算
    unique_x = unique(x_data);
    means = arrayfun(@(x) mean(y_data(x_data == x), 'omitnan'), unique_x);
    stds = arrayfun(@(x) std(y_data(x_data == x), 'omitnan'), unique_x);
    
    % 在第二个y轴上绘制平均长宽比
    yyaxis right;
    a = errorbar(unique_x, means, stds, 'ko-', 'LineWidth', 3, ...
                 'MarkerSize', 10, 'MarkerFaceColor', 'black', ...
                 "DisplayName", "Average value");
    ylabel('Average Aspect Ratio');
    set(gca, 'YColor', 'k');
    legend(a, 'Location', 'best');
    yyaxis left;
end

function save_result_plot(title)
    % 保存汇总统计图
    
    global ResultPath
    
    save_path = ResultPath;
    filename = title + ".jpg";
    full_path = fullfile(save_path, filename);
    print(gcf, full_path, '-djpeg', '-r300');
end

function plot_angle_distributions()
    % 绘制角度分布图
    
    global allRotatedPoints allAngles scan_parameter_values
    
    % 绘制旋转点的角度分布
    plotSeparatedAngleDistributions(allRotatedPoints, scan_parameter_values, 36);
    save_result_plot("NewtonPlot")
    
    
    % 绘制边角度的直方图
    plotAngleHistograms(allAngles, scan_parameter_values, 20);
    save_result_plot("AngleHistograms")
end