function visualizeComprehensiveBottleneckAnalysisWithSplitting(BW, CC, bottleneckResults, img, folderPath, fileNum)
% 可视化综合瓶颈分析结果
% 输入：
%   BW - 二值图像
%   CC - bwconncomp的输出
%   bottleneckResults - comprehensiveBottleneckAnalysis的输出结果
    
    numRegions = CC.NumObjects;

    cmap = zeros(numRegions + 1, 3); % +1 用于背景

    figure(12);
    imagesc(img, [0, 3])
    colormap jet
    clim([-0.5 3.5])
    hold on;

    labeled = labelmatrix(CC);
    % 获取所有区域的边界
    boundaries = bwboundaries(labeled, 8, 'noholes'); % 'noholes' 表示不包含孔洞
    
    for i = 1:numRegions
        hue = (i - 1) / max(1, numRegions - 1);
        cmap(i, :) = hsv2rgb([hue, 0.8, 0.9]);
        
        % 获取第i个区域的边界
        boundary = boundaries{i};
        
        % 注意：boundary返回的是[行, 列]，但plot需要的是[x, y] = [列, 行]
        cols = boundary(:, 2);
        rows = boundary(:, 1);
        
        % 绘制轮廓线
        plot(cols, rows, 'LineWidth', 2, 'Color', cmap(i,:));
    end
    
    hold off;

    % boundaries = bwboundaries(BW);
    % for k = 1:length(boundaries)
    %     boundary = boundaries{k};
    %     hue = (k - 1) / max(1, numRegions - 1);
    %     cmap(k, :) = hsv2rgb([hue, 0.8, 0.9]);
    %     plot(boundary(:,2), boundary(:,1), 'Color', cmap(k,:), 'LineWidth', 2.5);
    % end
    % 
    % legend('boundary');
    % hold off;
    
    % 子图3：显示距离变换和瓶颈位置
    % distTransform = bwdist(~BW);
    % imshow(distTransform, []);
    % hold on;
    
    % 标记有瓶颈的区域中心
    % for i = 1:numRegions
    %     if bottleneckResults(i).HasBottleneck
    %         % 获取区域属性
    %         regionBW = false(size(BW));
    %         regionBW(CC.PixelIdxList{i}) = true;
    %         stats = regionprops(regionBW, 'Centroid');
    % 
    %         if ~isempty(stats)
    %             centroid = stats.Centroid;
    %             plot(centroid(1), centroid(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    %         end
    %     end
    % end
    
    % labeled = labelmatrix(CC);
    % cmap = zeros(numRegions + 1, 3); % +1 用于背景
    % for i = 1:numRegions
    %     pixels = CC.PixelIdxList{i};
    %     [rows, cols] = ind2sub(size(img), pixels);
    %     hue = (i - 1) / max(1, numRegions - 1);
    %     cmap(i, :) = hsv2rgb([hue, 0.8, 0.9]);
    %     plot(cols, rows, '.', 'MarkerSize', 24, 'MarkerFaceColor', cmap(i,:));
    %     RGB = label2rgb(labeled, cmap, [0, 0, 0], 'shuffle');
    %     imagesc(RGB);
    % end

    % colorbar;
    % colormap gray

    hold off;

    save_path = folderPath + "AnalysisPlot\Split\";
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    filename = "Split" + fileNum + ".jpg";
    full_path = fullfile(save_path, filename);
    % 导出图像
    print(gcf, full_path, '-djpeg', '-r300'); 
end