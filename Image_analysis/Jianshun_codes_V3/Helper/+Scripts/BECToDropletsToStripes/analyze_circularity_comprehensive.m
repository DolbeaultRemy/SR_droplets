function [circularity_metrics, region_stats] = analyze_circularity_comprehensive(CC_split, selectFlag_row, show_plots)
    % 综合分析检测结构的圆度指标
    % 输入:
    %   CC_split - 连通域分割结构体
    %   selectFlag_row - 当前行的选择标志向量
    %   show_plots - 是否显示可视化图表 (可选，默认false)
    % 输出:
    %   circularity_metrics - 包含各种圆度指标的结构体
    %   region_stats - 区域统计信息
    
    % 设置默认参数
    if nargin < 3
        show_plots = false;
    end
    
    % 获取区域属性
    stats = regionprops(CC_split, 'Area', 'Perimeter', 'BoundingBox', ...
                        'MajorAxisLength', 'MinorAxisLength', ...
                        'Eccentricity', 'EquivDiameter', 'Solidity', ...
                        'Extent', 'Centroid');
    
    num_regions = length(stats);
    
    % 初始化输出结构
    circularity_metrics = struct();
    region_stats = stats;
    
    % 预分配数组
    area = zeros(num_regions, 1);
    perimeter = zeros(num_regions, 1);
    major_axis = zeros(num_regions, 1);
    minor_axis = zeros(num_regions, 1);
    eccentricity = zeros(num_regions, 1);
    equiv_diameter = zeros(num_regions, 1);
    solidity = zeros(num_regions, 1);
    extent = zeros(num_regions, 1);
    
    % 提取统计数据
    for i = 1:num_regions
        area(i) = stats(i).Area;
        perimeter(i) = stats(i).Perimeter;
        major_axis(i) = stats(i).MajorAxisLength;
        minor_axis(i) = stats(i).MinorAxisLength;
        eccentricity(i) = stats(i).Eccentricity;
        equiv_diameter(i) = stats(i).EquivDiameter;
        solidity(i) = stats(i).Solidity;
        extent(i) = stats(i).Extent;
    end
    
    % 计算各种圆度指标
    % 1. 基本圆度 (circularity)
    circularity = zeros(num_regions, 1);
    for i = 1:num_regions
        if perimeter(i) > 0 && area(i) > 0
            circularity(i) = (4 * pi * area(i)) / (perimeter(i)^2);
        else
            circularity(i) = 0;
        end
    end
    
    % 2. 圆整度 (roundness)
    roundness = zeros(num_regions, 1);
    for i = 1:num_regions
        if major_axis(i) > 0
            roundness(i) = (4 * area(i)) / (pi * major_axis(i)^2);
        else
            roundness(i) = 0;
        end
    end
    
    % 3. 紧密度 (compactness)
    compactness = zeros(num_regions, 1);
    for i = 1:num_regions
        if area(i) > 0
            compactness(i) = perimeter(i)^2 / area(i);
        else
            compactness(i) = Inf;
        end
    end
    
    % 4. 形状因子 (form factor)
    form_factor = zeros(num_regions, 1);
    for i = 1:num_regions
        if perimeter(i) > 0
            form_factor(i) = (4 * pi * area(i)) / (perimeter(i)^2);
        else
            form_factor(i) = 0;
        end
    end
    
    % 5. 纵横比 (aspect ratio)
    aspect_ratio = zeros(num_regions, 1);
    for i = 1:num_regions
        if minor_axis(i) > 0
            aspect_ratio(i) = major_axis(i) / minor_axis(i);
        else
            aspect_ratio(i) = 0;
        end
    end
    
    % 6. 椭圆度 (ellipticity)
    ellipticity = zeros(num_regions, 1);
    for i = 1:num_regions
        if minor_axis(i) > 0
            ellipticity(i) = major_axis(i) / minor_axis(i);
        else
            ellipticity(i) = 0;
        end
    end
    
    % 7. 圆度指数 (circularity index)
    circularity_index = zeros(num_regions, 1);
    for i = 1:num_regions
        if equiv_diameter(i) > 0 && perimeter(i) > 0
            circularity_index(i) = pi * equiv_diameter(i) / perimeter(i);
        else
            circularity_index(i) = 0;
        end
    end
    
    % 存储所有指标到输出结构
    circularity_metrics.area = area;
    circularity_metrics.perimeter = perimeter;
    circularity_metrics.major_axis = major_axis;
    circularity_metrics.minor_axis = minor_axis;
    circularity_metrics.eccentricity = eccentricity;
    circularity_metrics.equiv_diameter = equiv_diameter;
    circularity_metrics.solidity = solidity;
    circularity_metrics.extent = extent;
    circularity_metrics.circularity = circularity;
    circularity_metrics.roundness = roundness;
    circularity_metrics.compactness = compactness;
    circularity_metrics.form_factor = form_factor;
    circularity_metrics.aspect_ratio = aspect_ratio;
    circularity_metrics.ellipticity = ellipticity;
    circularity_metrics.circularity_index = circularity_index;
    
    % 根据选择标志筛选结果
    if nargin > 1 && ~isempty(selectFlag_row)
        selected_indices = find(selectFlag_row(1:num_regions));
        
        % 计算选中区域的统计摘要
        if ~isempty(selected_indices)
            fprintf('=== 圆度综合分析 ===\n');
            fprintf('选中区域数量: %d (总计: %d)\n', length(selected_indices), num_regions);
            fprintf('----------------------------------------\n');
            
            % 面积统计
            fprintf('面积统计:\n');
            fprintf('  平均面积: %.2f ± %.2f 像素\n', ...
                    mean(area(selected_indices)), std(area(selected_indices)));
            fprintf('  面积范围: [%.2f, %.2f] 像素\n', ...
                    min(area(selected_indices)), max(area(selected_indices)));
            
            % 主要圆度指标
            fprintf('\n主要圆度指标:\n');
            fprintf('  圆度(circularity): %.3f ± %.3f (理想圆形: 1.000)\n', ...
                    mean(circularity(selected_indices)), std(circularity(selected_indices)));
            fprintf('  圆整度(roundness): %.3f ± %.3f (理想圆形: 1.000)\n', ...
                    mean(roundness(selected_indices)), std(roundness(selected_indices)));
            fprintf('  紧密度(compactness): %.3f ± %.3f (理想圆形: %.3f)\n', ...
                    mean(compactness(selected_indices)), std(compactness(selected_indices)), 4*pi);
            
            % 形状特征
            fprintf('\n形状特征:\n');
            fprintf('  偏心率(eccentricity): %.3f ± %.3f (圆形: 0.000, 椭圆: 0.0-1.0)\n', ...
                    mean(eccentricity(selected_indices)), std(eccentricity(selected_indices)));
            fprintf('  纵横比(aspect ratio): %.3f ± %.3f (圆形: 1.000)\n', ...
                    mean(aspect_ratio(selected_indices)), std(aspect_ratio(selected_indices)));
            fprintf('  实心度(solidity): %.3f ± %.3f (完全实心: 1.000)\n', ...
                    mean(solidity(selected_indices)), std(solidity(selected_indices)));
            fprintf('  范围度(extent): %.3f ± %.3f (最大: 1.000)\n', ...
                    mean(extent(selected_indices)), std(extent(selected_indices)));
            
            % 分类统计
            classify_circularity_stats(circularity(selected_indices), roundness(selected_indices));
        end
    end
    
    % 可视化（如果需要）
    if show_plots && ~isempty(selected_indices)
        visualize_circularity_analysis(circularity_metrics, selected_indices);
    end
end

function classify_circularity_stats(circularity_values, roundness_values)
    % 根据圆度值对区域进行分类
    
    % 定义分类阈值
    thresholds = [0.9, 0.8, 0.7, 0.6, 0.5];
    labels = {'接近圆形', '较圆', '一般', '略不规则', '不规则'};
    
    fprintf('\n圆度分类统计:\n');
    fprintf('  圆度(circularity)分类:\n');
    
    for i = 1:length(thresholds)
        if i == 1
            count = sum(circularity_values >= thresholds(i));
            fprintf('    %s (≥%.2f): %d 个区域 (%.1f%%)\n', ...
                    labels{i}, thresholds(i), count, 100*count/length(circularity_values));
        else
            count = sum(circularity_values >= thresholds(i) & circularity_values < thresholds(i-1));
            fprintf('    %s (%.2f-%.2f): %d 个区域 (%.1f%%)\n', ...
                    labels{i}, thresholds(i), thresholds(i-1), count, 100*count/length(circularity_values));
        end
    end
    
    % 低于最低阈值的
    count = sum(circularity_values < thresholds(end));
    fprintf('    非常不规则 (<%.2f): %d 个区域 (%.1f%%)\n', ...
            thresholds(end), count, 100*count/length(circularity_values));
    
    % 圆形比例
    circular_count = sum(circularity_values >= 0.85);
    fprintf('  圆形比例 (圆度≥0.85): %.1f%%\n', 100*circular_count/length(circularity_values));
end

function visualize_circularity_analysis(circularity_metrics, selected_indices)
    % 可视化圆度分析结果
    
    % 创建图形窗口
    figure('Position', [100, 100, 1200, 800]);
    
    % 子图1: 圆度直方图
    subplot(2, 3, 1);
    histogram(circularity_metrics.circularity(selected_indices), 20, ...
              'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xlabel('圆度 (circularity)');
    ylabel('频率');
    title('圆度分布直方图');
    grid on;
    xlim([0, 1.2]);
    hold on;
    % 标记理想圆形位置
    plot([1, 1], ylim, 'r--', 'LineWidth', 2);
    text(1.05, max(ylim)*0.9, '完美圆形', 'Color', 'r', 'FontSize', 10);
    hold off;
    
    % 子图2: 圆度vs面积散点图
    subplot(2, 3, 2);
    scatter(circularity_metrics.area(selected_indices), ...
            circularity_metrics.circularity(selected_indices), 30, ...
            'filled', 'MarkerFaceAlpha', 0.6);
    xlabel('面积 (像素)');
    ylabel('圆度');
    title('圆度 vs 面积');
    grid on;
    hold on;
    % 添加趋势线
    x = circularity_metrics.area(selected_indices);
    y = circularity_metrics.circularity(selected_indices);
    valid_idx = ~isnan(x) & ~isnan(y);
    if sum(valid_idx) > 1
        p = polyfit(x(valid_idx), y(valid_idx), 1);
        x_fit = linspace(min(x(valid_idx)), max(x(valid_idx)), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    end
    hold off;
    
    % 子图3: 圆度vs偏心率散点图
    subplot(2, 3, 3);
    scatter(circularity_metrics.eccentricity(selected_indices), ...
            circularity_metrics.circularity(selected_indices), 30, ...
            'filled', 'MarkerFaceAlpha', 0.6);
    xlabel('偏心率 (eccentricity)');
    ylabel('圆度');
    title('圆度 vs 偏心率');
    grid on;
    hold on;
    % 理想圆形位置
    plot([0, 0], [0, 1], 'r--', 'LineWidth', 1);
    hold off;
    
    % 子图4: 不同圆度指标的箱线图
    subplot(2, 3, 4);
    data_to_plot = [circularity_metrics.circularity(selected_indices), ...
                    circularity_metrics.roundness(selected_indices)];
    boxplot(data_to_plot, 'Labels', {'圆度', '圆整度'});
    ylabel('指标值');
    title('圆度指标分布');
    grid on;
    
    % 子图5: 圆度vs纵横比散点图
    subplot(2, 3, 5);
    scatter(circularity_metrics.aspect_ratio(selected_indices), ...
            circularity_metrics.circularity(selected_indices), 30, ...
            'filled', 'MarkerFaceAlpha', 0.6);
    xlabel('纵横比 (aspect ratio)');
    ylabel('圆度');
    title('圆度 vs 纵横比');
    grid on;
    hold on;
    % 标记理想圆形位置
    plot([1, 1], [0, 1], 'r--', 'LineWidth', 1);
    hold off;
    
    % 子图6: 圆度分类饼图
    subplot(2, 3, 6);
    thresholds = [0.9, 0.8, 0.7, 0.6, 0.5];
    labels = {'接近圆形', '较圆', '一般', '略不规则', '不规则', '非常不规则'};
    
    circularity_values = circularity_metrics.circularity(selected_indices);
    counts = zeros(length(labels), 1);
    
    % 计算每类数量
    for i = 1:length(thresholds)
        if i == 1
            counts(i) = sum(circularity_values >= thresholds(i));
        else
            counts(i) = sum(circularity_values >= thresholds(i) & circularity_values < thresholds(i-1));
        end
    end
    counts(end) = sum(circularity_values < thresholds(end));
    
    % 绘制饼图
    pie(counts, labels);
    title('圆度分类比例');
    
    % 调整图形外观
    sgtitle('圆度综合分析可视化', 'FontSize', 16, 'FontWeight', 'bold');
    set(gcf, 'Color', 'w');
end