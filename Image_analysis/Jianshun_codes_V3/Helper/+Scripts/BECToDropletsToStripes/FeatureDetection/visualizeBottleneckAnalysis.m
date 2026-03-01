function visualizeBottleneckAnalysis(BW, patchProps, bottleneckWidths)
% 可视化瓶颈分析结果
% 输入：
%   BW - 二值图像
%   CC - bwconncomp的输出
%   bottleneckWidths - 方法一计算得到的瓶颈宽度数组

    figure;
    
    % 子图1：原始二值图像和区域标记
    subplot(2, 3, 1);
    imshow(BW);
    title('原始二值图像');
    
    % 子图2：彩色标记不同区域
    subplot(2, 3, 2);
    labeled = labelmatrix(patchProps);
    RGB = label2rgb(labeled, 'jet', 'k', 'shuffle');
    imshow(RGB);
    title('区域标记（不同颜色）');
    
    % 子图3：显示距离变换
    subplot(2, 3, 3);
    distTransform = bwdist(~BW);
    imshow(distTransform, []);
    colorbar;
    title('距离变换图');
    colormap(jet);
    
    % 子图4：显示骨架和瓶颈位置
    subplot(2, 3, 4);
    imshow(BW);
    hold on;
    
    % 为每个区域计算并绘制骨架和瓶颈点
    for i = 1:patchProps.NumObjects
        if bottleneckWidths(i) > 0
            % 创建当前区域的二值图像
            regionBW = false(size(BW));
            regionBW(patchProps.PixelIdxList{i}) = true;
            
            % 计算骨架
            skeleton = bwmorph(regionBW, 'skel', inf);
            
            % 计算距离变换
            regionDist = bwdist(~regionBW);
            
            % 找到骨架上的瓶颈点（最小宽度的点）
            skeletonWidths = 2 * regionDist(skeleton);
            minWidth = min(skeletonWidths);
            
            % 找到所有具有最小宽度的骨架点
            [skeletonY, skeletonX] = find(skeleton);
            skeletonPts = [skeletonX, skeletonY];
            widthsAtSkeleton = 2 * regionDist(skeleton);
            bottleneckIndices = find(widthsAtSkeleton == minWidth);
            
            if ~isempty(bottleneckIndices)
                % 绘制骨架
                plot(skeletonX, skeletonY, 'g.', 'MarkerSize', 3);
                
                % 绘制瓶颈点
                bottleneckX = skeletonX(bottleneckIndices);
                bottleneckY = skeletonY(bottleneckIndices);
                plot(bottleneckX, bottleneckY, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                
                % 标记瓶颈宽度
                for j = 1:length(bottleneckX)
                    text(bottleneckX(j)+5, bottleneckY(j)-5, ...
                        sprintf('%.1f', minWidth), ...
                        'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold');
                end
            end
        end
    end
    title('骨架（绿色）和瓶颈位置（红色圆圈）');
    
    % 子图5：显示瓶颈宽度分布
    subplot(2, 3, 5);
    validWidths = bottleneckWidths(bottleneckWidths > 0);
    if ~isempty(validWidths)
        histogram(validWidths, 20);
        xlabel('瓶颈宽度（像素）');
        ylabel('区域数量');
        title('瓶颈宽度分布');
        grid on;
    else
        text(0.5, 0.5, '无有效瓶颈数据', 'HorizontalAlignment', 'center');
        title('瓶颈宽度分布');
    end
    
    % 子图6：区域统计信息
    subplot(2, 3, 6);
    % 创建统计表格
    statsTable = cell(patchProps.NumObjects+1, 3);
    statsTable{1,1} = '区域编号';
    statsTable{1,2} = '瓶颈宽度';
    statsTable{1,3} = '是否有瓶颈';
    
    for i = 1:patchProps.NumObjects
        statsTable{i+1,1} = sprintf('%d', i);
        statsTable{i+1,2} = sprintf('%.2f', bottleneckWidths(i));
        if bottleneckWidths(i) > 0
            % 简单判断是否有瓶颈：如果宽度小于区域等效直径的1/3则认为有瓶颈
            regionArea = length(patchProps.PixelIdxList{i});
            equivalentDiameter = 2 * sqrt(regionArea/pi);
            if bottleneckWidths(i) < equivalentDiameter / 3
                statsTable{i+1,3} = '是';
            else
                statsTable{i+1,3} = '否';
            end
        else
            statsTable{i+1,3} = '否';
        end
    end
    
    % 显示表格
    axis off;
    text(0, 0.9, '区域瓶颈分析统计', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 创建表格文本
    tableY = 0.8;
    lineHeight = 0.06;
    for row = 1:size(statsTable,1)
        for col = 1:size(statsTable,2)
            textX = (col-1)*0.3;
            if row == 1
                % 表头加粗
                text(textX, tableY, statsTable{row,col}, ...
                    'FontWeight', 'bold', 'FontSize', 10);
            else
                % 数据行
                if col == 3 && strcmp(statsTable{row,col}, '是')
                    % 有瓶颈的标记为红色
                    text(textX, tableY, statsTable{row,col}, ...
                        'Color', 'red', 'FontSize', 9);
                else
                    text(textX, tableY, statsTable{row,col}, 'FontSize', 9);
                end
            end
        end
        tableY = tableY - lineHeight;
    end
    
    % 添加整体统计信息
    numBottlenecks = sum(cellfun(@(x) strcmp(x, '是'), statsTable(2:end,3)));
    text(0, tableY - 0.1, sprintf('总结: %d/%d 个区域存在瓶颈结构', ...
        numBottlenecks, patchProps.NumObjects), 'FontSize', 11, 'FontWeight', 'bold');
    
    % 设置图形属性
    set(gcf, 'Position', [100, 100, 1400, 800]);
    sgtitle('区域瓶颈结构分析可视化结果', 'FontSize', 14, 'FontWeight', 'bold');
end

% 使用示例
function runBottleneckAnalysisExample()
    % 创建示例图像或读取您的图像
    % 这里创建一个包含不同形状的示例图像
    BW = false(200, 200);
    
    % 添加一些测试形状
    % 圆形（无瓶颈）
    [x, y] = meshgrid(1:200);
    circle = (x-50).^2 + (y-50).^2 <= 400;
    BW = BW | circle;
    
    % 哑铃形状（有瓶颈）
    BW(80:120, 80:85) = true;  % 左端
    BW(80:120, 115:120) = true; % 右端
    BW(95:105, 85:115) = true;  % 狭窄部分
    
    % 不规则形状
    BW(150:180, 30:60) = true;
    BW(160:170, 60:80) = true;
    BW(150:180, 80:100) = true;
    
    % 连通组件分析
    CC = bwconncomp(BW);
    
    % 使用方法一计算瓶颈宽度
    bottleneckWidths = findBottleneckWidths(BW, CC);
    
    % 可视化结果
    visualizeBottleneckAnalysis(BW, CC, bottleneckWidths);
end

% 您的主程序中使用：
% 假设您已经有了 BW 和 CC
% bottleneckWidths = findBottleneckWidths(BW, CC);
% visualizeBottleneckAnalysis(BW, CC, bottleneckWidths);