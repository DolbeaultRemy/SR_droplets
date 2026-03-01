function visualizeComprehensiveBottleneckAnalysisWithSplitting(BW, CC, bottleneckResults, img, folderPath, fileName)
% VISUALIZECOMPREHENSIVEBOTTLENECKANALYSISWITHSPLITTING
%   Display the split regions with distinct colors.
% Inputs:
%   BW                - binary image (not used directly, but kept for consistency)
%   CC                - bwconncomp after splitting
%   bottleneckResults - struct array from analysis (not used in this simple version)
%   img               - original cropped image for background
%   folderPath        - string to save figure
%   fileName          - base name for saved image

    numRegions = CC.NumObjects;
    labeled = labelmatrix(CC);
    % Create a colormap with distinct colors for each region
    cmap = zeros(numRegions, 3);
    for i = 1:numRegions
        hue = (i - 1) / max(1, numRegions - 1);
        cmap(i, :) = hsv2rgb([hue, 0.8, 0.9]);
    end
    % Convert labeled image to RGB
    rgbImage = label2rgb(labeled, cmap, [0,0,0]);

    figure(12); clf;
    % Show original image in grayscale as background
    imagesc(img, [0, 3]);
    colormap(gca, 'gray');
    hold on;
    % Overlay the colored regions with some transparency
    h = imshow(rgbImage);
    set(h, 'AlphaData', 0.5); % transparency
    % Plot boundaries for clarity
    boundaries = bwboundaries(labeled, 8, 'noholes');
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'Color', cmap(k,:), 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Split into %d Regions', numRegions));

    % Save figure
    savePath = fullfile(folderPath, 'AnalysisPlot', 'Split');
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    saveas(gcf, fullfile(savePath, ['Split_' fileName '.jpg']));
    % Alternatively use print for higher resolution:
    % print(gcf, fullfile(savePath, ['Split_' fileName '.jpg']), '-djpeg', '-r300');
end