% plotSeparatedAngleDistributions.m
function plotSeparatedAngleDistributions(allRotatedPoints, scan_parameter_values, numBins)
    % 将每个扫描参数的角度分布绘制在不同的子图中
    % 输入: allRotatedPoints - 三维数组 (scan_params × points × 2)
    %       scan_parameter_values - 扫描参数值数组
    %       numBins - 角度区间的数量
    
    if nargin < 3
        numBins = 36; % 默认36个扇区，每10度一个
    end
    
    % 获取唯一扫描参数值
    [unique_scans, ~, idx_to_unique] = unique(scan_parameter_values);
    n_unique_scans = length(unique_scans);
    
    fprintf('发现 %d 个唯一扫描参数（原始 %d 个参数）\n', ...
            n_unique_scans, length(scan_parameter_values));
    
    % 计算子图布局
    [nRows, nCols] = calculate_subplot_layout(n_unique_scans);
    
    % 创建图形
    figure('Position', [100, 100, 1600, 900]);
    
    % 生成颜色映射
    colors = generate_colormap(n_unique_scans);
    
    % 为每个唯一参数创建子图
    for i = 1:n_unique_scans
        create_single_angle_subplot(i, unique_scans, idx_to_unique, ...
                                    allRotatedPoints, numBins, colors, ...
                                    nRows, nCols);
    end
    
    theme("light");
    fprintf('\n绘图完成: %d 个子图已创建\n', n_unique_scans);
end

function [nRows, nCols] = calculate_subplot_layout(n_unique_scans)
    % 计算子图布局
    
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
end

function colors = generate_colormap(n_unique_scans)
    % 生成颜色映射
    
    if n_unique_scans <= 10
        % 对于少数参数，使用预定义的鲜艳颜色
        colors = lines(n_unique_scans);
    else
        % 对于多个参数，使用jet颜色映射
        colors = jet(n_unique_scans);
    end
end

function create_single_angle_subplot(idx, unique_scans, idx_to_unique, ...
                                     allRotatedPoints, numBins, colors, ...
                                     nRows, nCols)
    % 创建单个角度分布子图
    
    % 选择当前子图
    ax = subplot(nRows, nCols, idx, polaraxes);
    
    % 找到所有等于当前唯一值的索引
    same_param_idx = find(idx_to_unique == idx);
    
    % 合并这些索引对应的所有点
    temp_points = [];
    for j = 1:length(same_param_idx)
        file_idx = same_param_idx(j);
        current_points = squeeze(allRotatedPoints(file_idx, :, :));
        
        % 移除无效数据
        valid_idx = all(isfinite(current_points), 2);
        valid_idx(all(current_points == 0, 2)) = false;
        current_points = current_points(valid_idx, :);
        
        % 添加到临时存储
        temp_points = [temp_points; current_points];
    end
    
    % 如果有点数据，绘制分布
    if ~isempty(temp_points)
        plot_angle_distribution(ax, temp_points, numBins, colors(idx, :), ...
                                unique_scans(idx));
    else
        % 创建空极坐标图
        title(sprintf('B = %.4f (G)\nNo Data', unique_scans(idx)), ...
              'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
        text(0.5, 0.5, 'No Data', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
    end
    
    thetalim([0, 180]);  % 只显示0°到180°
    set(gca, 'LineWidth', 1.5, 'Box', 'on');
end

function plot_angle_distribution(ax, points, numBins, color, scan_value)
    % 绘制角度分布
    
    % 分离x和y坐标
    x_data = points(:, 1);
    y_data = points(:, 2);
    
    % 计算每个点的角度
    angles_rad = atan2(y_data, x_data);
    angles_deg = mod(rad2deg(angles_rad), 360);
    
    % 计算角度直方图
    angle_edges = linspace(0, 360, numBins + 1);
    angle_centers = angle_edges(1:end-1) + (angle_edges(2) - angle_edges(1))/2;
    counts = histcounts(angles_deg, angle_edges);
    
    % 处理对称性
    counts = counts + fliplr(counts);
    
    % 将角度转换为弧度
    theta_rad = deg2rad(angle_centers);
    theta_plot = theta_rad;
    r_plot = counts;
    
    % 创建极坐标图
    polarplot(ax, theta_plot, r_plot, '-', 'Color', color, 'LineWidth', 2);
    hold on;
    
    % 设置极坐标图属性
    ax.ThetaZeroLocation = 'right';
    ax.ThetaDir = 'counterclockwise';
    ax.RGrid = 'on';
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.3;
    ax.ThetaTick = 0:30:330;
    
    % 添加标题
    total_points = length(angles_deg);
    title_str = sprintf('B = %.4f(G), n = %d', scan_value, total_points);
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    fprintf('参数 %d (值=%.4f): %d 个点\n', ...
            scan_value, total_points);
end