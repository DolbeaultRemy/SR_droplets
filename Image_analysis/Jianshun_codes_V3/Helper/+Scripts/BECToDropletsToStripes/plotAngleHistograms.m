% plotAngleHistograms.m
function plotAngleHistograms(allAngles, scan_parameter_values, numBins)
    % 将每个扫描参数的角度直方图绘制在不同的子图中
    % 输入: allAngles - 二维数组 (scan_params × data_points)
    %       scan_parameter_values - 扫描参数值数组
    %       numBins - 直方图的箱数
    
    if nargin < 3
        numBins = 18; % 默认18个箱，每10度一个
    end
    
    % 验证输入数据维度
    n_scans = length(scan_parameter_values);
    if size(allAngles, 1) ~= n_scans
        error('allAngles的第一个维度大小(%d)与scan_parameter_values长度(%d)不匹配', ...
              size(allAngles, 1), n_scans);
    end
    
    % 获取数组维度信息
    [n_scans, n_points] = size(allAngles);
    fprintf('数组维度: %d × %d\n', n_scans, n_points);
    
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
        create_single_histogram_subplot(i, unique_scans, idx_to_unique, ...
                                        allAngles, numBins, colors, ...
                                        nRows, nCols);
    end
    
    % 为整个图形添加标题
    sgtitle('Angular distribution histogram', 'FontSize', 14, 'FontWeight', 'bold');
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

function create_single_histogram_subplot(idx, unique_scans, idx_to_unique, ...
                                         allAngles, numBins, colors, ...
                                         nRows, nCols)
    % 创建单个直方图子图
    
    ax = subplot(nRows, nCols, idx);
    
    % 找到所有等于当前唯一值的索引
    same_param_idx = find(idx_to_unique == idx);
    
    % 合并这些索引对应的所有角度数据
    temp_angles = [];
    for j = 1:length(same_param_idx)
        file_idx = same_param_idx(j);
        current_angles = allAngles(file_idx, :);
        
        % 移除无效数据
        valid_idx = isfinite(current_angles) & ...
                   (current_angles > 0) & ...     % 排除0值
                   (current_angles <= 360);       % 只保留0-360度的值
        valid_angles = current_angles(valid_idx);
        
        % 添加到总的角度数据
        temp_angles = [temp_angles, valid_angles];
    end
    
    % 如果有角度数据，绘制直方图
    if ~isempty(temp_angles)
        plot_angle_histogram(ax, temp_angles, numBins, colors(idx, :), ...
                             unique_scans(idx));
    else
        % 创建空图
        title(sprintf('α = %.4f\nNo Data', unique_scans(idx)), ...
              'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
        text(0.5, 0.5, 'No Data', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
        xlim([0, 180]);
        ylim([0, 1]);
        xlabel('Degree', 'FontSize', 10);
        ylabel('count', 'FontSize', 10);
        xticks(0:30:180);
    end
    
    set(gca, 'LineWidth', 1, 'Box', 'on');
end

function plot_angle_histogram(ax, angles, numBins, color, scan_value)
    % 绘制角度直方图
    
    % 将角度归一化到0-360度范围
    angles_deg = mod(angles, 360);
    angles_deg(angles_deg == 0) = 360;  % 处理边界情况
    
    % 映射到0-180度范围（对称处理）
    mapped_angles = angles_deg;
    mapped_angles(mapped_angles > 180) = 360 - mapped_angles(mapped_angles > 180);
    
    % 计算直方图
    [counts, bin_edges] = histcounts(mapped_angles, numBins, ...
                                     'BinEdges', 0:180/numBins:180);
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    % 归一化计数
    counts = counts / max(counts);
    
    % 绘制条形图
    bar(ax, bin_centers, counts, 'FaceColor', color, 'EdgeColor', color, ...
        'FaceAlpha', 0.7);
    hold on;
    
    % 添加平滑曲线
    smooth_x = linspace(0, 180, 500);
    smooth_counts = interp1([bin_centers(1)-180/numBins, bin_centers, ...
                             bin_centers(end)+180/numBins], ...
                           [counts(end), counts, counts(1)], smooth_x, 'pchip');
    plot(ax, smooth_x, smooth_counts, 'Color', color, 'LineWidth', 2);
    
    % 设置坐标轴属性
    xlim([0, 180]);
    xlabel('Degree', 'FontSize', 10);
    ylabel('Normalized Count', 'FontSize', 10);
    xticks(0:30:180);
    xticklabels({'0°', '30°', '60°', '90°', '120°', '150°', '180°'});
    
    % 添加网格
    grid on;
    grid minor;
    
    % 添加标题和信息
    title_str = sprintf('B = %.4f (G), n = %d', scan_value, length(mapped_angles));
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    fprintf('参数 %d (值=%.4f): %d 个角度值 (排除0值后)\n', ...
            scan_value, length(mapped_angles));
end