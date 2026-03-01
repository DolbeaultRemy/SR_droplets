function visualizeComprehensiveBottleneckAnalysis(BW, CC, bottleneckResults)
% 可视化综合瓶颈分析结果
% 输入：
%   BW - 二值图像
%   CC - bwconncomp的输出
%   bottleneckResults - comprehensiveBottleneckAnalysis的输出结果

    figure;
    
    % 子图1：原始二值图像和区域标记
    subplot(2, 3, 1);
    imshow(BW);
    title('原始二值图像');
    
    % 子图2：彩色标记不同区域，用红色标记有瓶颈的区域
    subplot(2, 3, 2);
    labeled = labelmatrix(CC);
    
    % 创建自定义颜色映射：有瓶颈的区域用红色，无瓶颈的用其他颜色
    hasBottleneck = [bottleneckResults.HasBottleneck];
    numRegions = CC.NumObjects;
    
    % 生成颜色映射
    cmap = zeros(numRegions + 1, 3); % +1 用于背景
    cmap(1, :) = [0, 0, 0]; % 背景黑色
    
    for i = 1:numRegions
        if hasBottleneck(i)
            cmap(i + 1, :) = [1, 0, 0]; % 有瓶颈：红色
        else
            % 无瓶颈：使用jet色彩空间的不同颜色
            hue = (i - 1) / max(1, numRegions - 1);
            cmap(i + 1, :) = hsv2rgb([hue, 0.8, 0.9]);
        end
    end
    
    RGB = label2rgb(labeled, cmap, [0, 0, 0], 'shuffle');
    imshow(RGB);
    title('区域标记（红色=有瓶颈）');
    
    % 子图3：显示距离变换和瓶颈位置
    subplot(2, 3, 3);
    distTransform = bwdist(~BW);
    imshow(distTransform, []);
    hold on;
    
    % 标记有瓶颈的区域中心
    for i = 1:numRegions
        if bottleneckResults(i).HasBottleneck
            % 获取区域属性
            regionBW = false(size(BW));
            regionBW(CC.PixelIdxList{i}) = true;
            stats = regionprops(regionBW, 'Centroid');
            
            if ~isempty(stats)
                centroid = stats.Centroid;
                plot(centroid(1), centroid(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
                text(centroid(1)+5, centroid(2)-5, ...
                    sprintf('%.1f', bottleneckResults(i).MinWidth), ...
                    'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold', ...
                    'BackgroundColor', 'white');
            end
        end
    end
    colorbar;
    title('距离变换 + 瓶颈标记');
    colormap(jet);
    
    % 子图4：显示骨架和详细瓶颈信息
    subplot(2, 3, 4);
    imshow(BW);
    hold on;
    
    % 为每个有瓶颈的区域绘制详细分析
    for i = 1:numRegions
        if bottleneckResults(i).HasBottleneck
            % 创建当前区域的二值图像
            regionBW = false(size(BW));
            regionBW(CC.PixelIdxList{i}) = true;
            
            % 计算骨架
            skeleton = bwmorph(regionBW, 'skel', inf);
            
            % 计算距离变换
            regionDist = bwdist(~regionBW);
            
            % 找到骨架上的瓶颈点（最小宽度的点）
            [skeletonY, skeletonX] = find(skeleton);
            widthsAtSkeleton = 2 * regionDist(skeleton);
            minWidth = bottleneckResults(i).MinWidth;
            bottleneckIndices = find(abs(widthsAtSkeleton - minWidth) < 0.1); % 容差
            
            if ~isempty(bottleneckIndices)
                % 绘制骨架
                plot(skeletonX, skeletonY, 'g.', 'MarkerSize', 3);
                
                % 绘制瓶颈点
                bottleneckX = skeletonX(bottleneckIndices);
                bottleneckY = skeletonY(bottleneckIndices);
                plot(bottleneckX, bottleneckY, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                
                % 标记瓶颈详细信息
                if ~isempty(bottleneckX)
                    text(bottleneckX(1)+5, bottleneckY(1)-5, ...
                        sprintf('min=%.1f\nmax=%.1f\nratio=%.2f', ...
                        bottleneckResults(i).MinWidth, ...
                        bottleneckResults(i).MaxWidth, ...
                        bottleneckResults(i).BottleneckRatio), ...
                        'Color', 'red', 'FontSize', 8, 'FontWeight', 'bold', ...
                        'BackgroundColor', 'white', 'VerticalAlignment', 'top');
                end
            end
        end
    end
    title('骨架和瓶颈详细信息');
    
    % 子图5：显示瓶颈比例分布
    subplot(2, 3, 5);
    bottleneckRatios = [bottleneckResults.BottleneckRatio];
    hasBottleneckIdx = find([bottleneckResults.HasBottleneck]);
    
    if ~isempty(hasBottleneckIdx)
        % 绘制有瓶颈区域的宽度比例
        histogram(bottleneckRatios(hasBottleneckIdx), 20, 'FaceColor', 'red');
        hold on;
        
        % 绘制所有区域的宽度比例
        histogram(bottleneckRatios, 20, 'FaceColor', 'blue', 'FaceAlpha', 0.3);
        
        xlabel('瓶颈比例 (MinWidth/MaxWidth)');
        ylabel('区域数量');
        title('瓶颈比例分布');
        legend('有瓶颈区域', '所有区域', 'Location', 'best');
        grid on;
    else
        text(0.5, 0.5, '无瓶颈区域', 'HorizontalAlignment', 'center');
        title('瓶颈比例分布');
    end
    
    % 子图6：详细统计信息表格
    subplot(2, 3, 6);
    axis off;
    
    % 创建详细的统计表格
    statsTable = cell(numRegions + 1, 6);
    statsTable{1,1} = '区域';
    statsTable{1,2} = '瓶颈';
    statsTable{1,3} = '最小宽';
    statsTable{1,4} = '最大宽';
    statsTable{1,5} = '比例';
    statsTable{1,6} = '面积';
    
    for i = 1:numRegions
        statsTable{i+1,1} = sprintf('%d', i);
        statsTable{i+1,2} = boolToText(bottleneckResults(i).HasBottleneck);
        statsTable{i+1,3} = sprintf('%.2f', bottleneckResults(i).MinWidth);
        statsTable{i+1,4} = sprintf('%.2f', bottleneckResults(i).MaxWidth);
        statsTable{i+1,5} = sprintf('%.3f', bottleneckResults(i).BottleneckRatio);
        statsTable{i+1,6} = sprintf('%d', bottleneckResults(i).Area);
    end
    
    % 显示表格
    text(0, 0.95, '详细区域统计', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 创建表格文本
    tableY = 0.85;
    lineHeight = 0.05;
    colWidths = [0.1, 0.15, 0.2, 0.2, 0.2, 0.15]; % 各列宽度
    
    for row = 1:size(statsTable,1)
        xPos = 0;
        for col = 1:size(statsTable,2)
            if row == 1
                % 表头加粗
                text(xPos, tableY, statsTable{row,col}, ...
                    'FontWeight', 'bold', 'FontSize', 9);
            else
                % 数据行
                textColor = 'black';
                if col == 2 && strcmp(statsTable{row,col}, '是')
                    textColor = 'red';
                elseif col == 5
                    % 根据瓶颈比例着色
                    ratio = bottleneckResults(row-1).BottleneckRatio;
                    if ratio < 0.2
                        textColor = [0.8, 0, 0]; % 深红
                    elseif ratio < 0.4
                        textColor = [1, 0.4, 0]; % 橙色
                    end
                end
                
                text(xPos, tableY, statsTable{row,col}, ...
                    'FontSize', 8, 'Color', textColor);
            end
            xPos = xPos + colWidths(col);
        end
        tableY = tableY - lineHeight;
    end
    
    % 添加整体统计信息
    numBottlenecks = sum([bottleneckResults.HasBottleneck]);
    avgRatio = mean([bottleneckResults.BottleneckRatio]);
    minRatio = min([bottleneckResults.BottleneckRatio]);
    
    summaryText = {...
        sprintf('区域总数: %d', numRegions), ...
        sprintf('有瓶颈区域: %d (%.1f%%)', numBottlenecks, numBottlenecks/numRegions*100), ...
        sprintf('平均瓶颈比例: %.3f', avgRatio), ...
        sprintf('最小瓶颈比例: %.3f', minRatio) ...
    };
    
    for i = 1:length(summaryText)
        text(0, tableY - i*0.04, summaryText{i}, 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % 设置图形属性
    set(gcf, 'Position', [100, 100, 1400, 800]);
    sgtitle('综合瓶颈分析可视化结果', 'FontSize', 14, 'FontWeight', 'bold');
end

function text = boolToText(boolValue)
% 将布尔值转换为中文文本
    if boolValue
        text = '是';
    else
        text = '否';
    end
end

% 可选：单独显示每个有瓶颈区域的详细分析
function visualizeIndividualBottleneckRegions(BW, CC, bottleneckResults)
% 为每个有瓶颈的区域单独显示详细分析
    
    hasBottleneck = [bottleneckResults.HasBottleneck];
    bottleneckIndices = find(hasBottleneck);
    
    for idx = 1:length(bottleneckIndices)
        i = bottleneckIndices(idx);
        
        figure;
        
        % 创建当前区域的二值图像
        regionBW = false(size(BW));
        regionBW(CC.PixelIdxList{i}) = true;
        
        % 子图1：区域二值图像和基本信息
        subplot(2, 3, 1);
        imshow(regionBW);
        title(sprintf('区域 %d', i));
        
        % 添加基本信息文本
        infoText = {...
            sprintf('面积: %d 像素', bottleneckResults(i).Area), ...
            sprintf('最小宽度: %.2f', bottleneckResults(i).MinWidth), ...
            sprintf('最大宽度: %.2f', bottleneckResults(i).MaxWidth), ...
            sprintf('平均宽度: %.2f', bottleneckResults(i).MeanWidth), ...
            sprintf('瓶颈比例: %.3f', bottleneckResults(i).BottleneckRatio) ...
        };
        
        for j = 1:length(infoText)
            text(10, 15 + j*20, infoText{j}, 'Color', 'white', ...
                'FontSize', 10, 'FontWeight', 'bold', ...
                'BackgroundColor', 'black');
        end
        
        % 子图2：距离变换热图
        subplot(2, 3, 2);
        regionDist = bwdist(~regionBW);
        imshow(regionDist, []);
        colorbar;
        title('距离变换');
        colormap(jet);
        
        % 子图3：骨架和瓶颈点
        subplot(2, 3, 3);
        imshow(regionBW);
        hold on;
        
        % 计算骨架
        skeleton = bwmorph(regionBW, 'skel', inf);
        [skeletonY, skeletonX] = find(skeleton);
        
        % 找到瓶颈点
        widthsAtSkeleton = 2 * regionDist(skeleton);
        minWidth = bottleneckResults(i).MinWidth;
        bottleneckIndices = find(abs(widthsAtSkeleton - minWidth) < 0.1);
        
        % 绘制骨架和瓶颈点
        plot(skeletonX, skeletonY, 'g.', 'MarkerSize', 5);
        if ~isempty(bottleneckIndices)
            bottleneckX = skeletonX(bottleneckIndices);
            bottleneckY = skeletonY(bottleneckIndices);
            plot(bottleneckX, bottleneckY, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            
            % 标记瓶颈宽度
            for j = 1:min(3, length(bottleneckX)) % 最多标记3个点
                text(bottleneckX(j)+3, bottleneckY(j)-3, ...
                    sprintf('%.1f', minWidth), ...
                    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'BackgroundColor', 'white');
            end
        end
        title('骨架和瓶颈点');
        
        % 子图4：宽度分布直方图
        subplot(2, 3, 4);
        histogram(widthsAtSkeleton, 30, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
        hold on;
        
        % 标记关键统计量
        ylims = ylim;
        plot([minWidth, minWidth], ylims, 'r-', 'LineWidth', 2);
        plot([bottleneckResults(i).MaxWidth, bottleneckResults(i).MaxWidth], ylims, 'g-', 'LineWidth', 2);
        plot([bottleneckResults(i).MeanWidth, bottleneckResults(i).MeanWidth], ylims, 'm-', 'LineWidth', 2);
        
        xlabel('宽度（像素）');
        ylabel('骨架点数');
        title('宽度分布');
        legend('宽度分布', '最小宽度', '最大宽度', '平均宽度', 'Location', 'best');
        grid on;
        
        % 子图5：区域边界和Feret直径示意
        subplot(2, 3, 5);
        imshow(regionBW);
        hold on;
        
        % 绘制边界
        boundaries = bwboundaries(regionBW);
        if ~isempty(boundaries)
            boundary = boundaries{1};
            plot(boundary(:,2), boundary(:,1), 'y-', 'LineWidth', 1);
        end
        
        % 标记最小宽度位置（近似）
        if ~isempty(bottleneckX)
            % 在第一个瓶颈点处绘制一条垂直线表示宽度
            x = bottleneckX(1);
            y = bottleneckY(1);
            plot([x, x], [y-5, y+5], 'r-', 'LineWidth', 3);
            text(x+5, y, sprintf('≈%.1f', minWidth*2), ...
                'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold');
        end
        title('边界和宽度示意');
        
        % 子图6：瓶颈严重程度评估
        subplot(2, 3, 6);
        axis off;
        
        ratio = bottleneckResults(i).BottleneckRatio;
        if ratio < 0.1
            severity = '严重瓶颈';
            color = [0.8, 0, 0];
            description = '区域存在明显的狭窄结构';
        elseif ratio < 0.3
            severity = '中度瓶颈';
            color = [1, 0.5, 0];
            description = '区域存在可识别的狭窄';
        elseif ratio < 0.5
            severity = '轻度瓶颈';
            color = [1, 0.8, 0];
            description = '区域宽度变化较明显';
        else
            severity = '无显著瓶颈';
            color = [0, 0.6, 0];
            description = '区域宽度相对均匀';
        end
        
        text(0.1, 0.9, '瓶颈严重程度评估', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.1, 0.7, sprintf('等级: %s', severity), 'FontSize', 14, ...
            'FontWeight', 'bold', 'Color', color);
        text(0.1, 0.5, description, 'FontSize', 10, 'VerticalAlignment', 'top');
        text(0.1, 0.3, sprintf('瓶颈比例: %.3f', ratio), 'FontSize', 11);
        
        set(gcf, 'Position', [100, 100, 1200, 700]);
        sgtitle(sprintf('区域 %d 详细瓶颈分析', i), 'FontSize', 16, 'FontWeight', 'bold');
    end
end

% 使用示例
function runComprehensiveExample()
    % 创建测试数据
    BW = false(200, 200);
    
    % 添加各种测试形状
    % 圆形（无瓶颈）
    [x, y] = meshgrid(1:200);
    circle1 = (x-50).^2 + (y-50).^2 <= 400;
    circle2 = (x-150).^2 + (y-50).^2 <= 250;
    BW = BW | circle1 | circle2;
    
    % 哑铃形状（有瓶颈）
    BW(80:120, 80:85) = true;
    BW(80:120, 115:120) = true;
    BW(95:105, 85:115) = true;
    
    % 不规则形状（可能有瓶颈）
    BW(140:170, 140:160) = true;
    BW(150:160, 160:180) = true;
    BW(140:170, 180:200) = true;
    
    % 连通组件分析
    CC = bwconncomp(BW);
    
    % 进行综合瓶颈分析
    bottleneckResults = comprehensiveBottleneckAnalysis(BW, CC, 0.3);
    
    % 可视化结果
    visualizeComprehensiveBottleneckAnalysis(BW, CC, bottleneckResults);
    
    % 可选：单独显示有瓶颈的区域
    % visualizeIndividualBottleneckRegions(BW, CC, bottleneckResults);
end