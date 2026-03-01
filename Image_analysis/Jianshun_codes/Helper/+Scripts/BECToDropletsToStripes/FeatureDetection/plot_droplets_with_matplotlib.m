function fig = plot_droplets_with_matplotlib(image_array, droplets, varargin)
% 在图像上绘制检测到的液滴
% 参数：
%     image_array: 二维或三维图像数组
%     droplets: 液滴列表，每个元素为[x, y, radius]向量
%     可选参数：
%         edge_color: 边缘颜色（默认'red'）
%         center_color: 中心颜色（默认'red'）
%         edge_linewidth: 边缘线宽（默认2）
%         center_size: 中心点大小（默认20）
%         figsize: 图像尺寸[宽,高]（默认[10,8]，单位英寸）
%         dpi: 输出分辨率（默认100）
% 返回：
%     fig: 图形句柄

% 解析输入参数
p = inputParser;
addRequired(p, 'image_array', @(x) isnumeric(x));
addRequired(p, 'droplets', @(x) iscell(x) || (isnumeric(x) && size(x,2)==3));
addParameter(p, 'edge_color', 'red', @(x) ischar(x) || (isnumeric(x) && length(x)==3));
addParameter(p, 'center_color', 'red', @(x) ischar(x) || (isnumeric(x) && length(x)==3));
addParameter(p, 'edge_linewidth', 2, @isnumeric);
addParameter(p, 'center_size', 20, @isnumeric);
addParameter(p, 'figsize', [10, 8], @(x) isnumeric(x) && length(x)==2);
addParameter(p, 'dpi', 100, @isnumeric);

parse(p, image_array, droplets, varargin{:});

% 创建图形
fig = figure('Units', 'inches', 'Position', [1, 1, p.Results.figsize]);
ax = axes(fig);

% 设置DPI（在MATLAB中保存时使用）
set(fig, 'PaperPositionMode', 'auto');

% 显示图像
if ndims(image_array) == 2
    % 灰度图像
    im = imshow(image_array, 'Parent', ax);
    colormap(ax, 'jet');
    caxis(ax, [0, 1]);
else
    % 彩色图像
    if size(image_array, 3) == 3
        % 如果是BGR格式需要转换为RGB
        rgb_image = image_array(:, :, [3, 2, 1]); % BGR转RGB
        im = imshow(rgb_image, 'Parent', ax);
    else
        im = imshow(image_array, 'Parent', ax);
    end
end

hold(ax, 'on');

% 确保droplets是cell数组格式
if isnumeric(droplets)
    temp_droplets = cell(size(droplets, 1), 1);
    for i = 1:size(droplets, 1)
        temp_droplets{i} = droplets(i, :);
    end
    droplets = temp_droplets;
end

% 绘制每个液滴
for i = 1:length(droplets)
    droplet = droplets{i};
    x = droplet(1);
    y = droplet(2);
    r = droplet(3);
    
    % 添加边缘圆环
    theta = 0:0.01:2*pi;
    x_circle = x + r * cos(theta);
    y_circle = y + r * sin(theta);
    
    plot(ax, x_circle, y_circle, ...
         'Color', p.Results.edge_color, ...
         'LineWidth', p.Results.edge_linewidth, ...
         'LineStyle', '-');
    
    % 添加中心点
    if strcmp(p.Results.center_color, 'white')
        scatter(ax, x, y, p.Results.center_size, ...
                'MarkerFaceColor', p.Results.center_color, ...
                'MarkerEdgeColor', 'black', ...
                'LineWidth', 1);
    else
        scatter(ax, x, y, p.Results.center_size, ...
                'MarkerFaceColor', p.Results.center_color, ...
                'MarkerEdgeColor', p.Results.center_color);
    end
end

% 设置坐标轴属性
axis(ax, 'image'); % 保持纵横比
set(ax, 'XTick', [], 'YTick', []);
set(ax, 'Color', 'none'); % 透明背景

% 自动调整布局
set(ax, 'Position', [0, 0, 1, 1]); % 使坐标轴填满图形

hold(ax, 'off');

end