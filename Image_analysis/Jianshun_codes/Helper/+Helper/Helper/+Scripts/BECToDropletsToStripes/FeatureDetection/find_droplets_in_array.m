function droplets = find_droplets_in_array(image_array, density_thresh, min_radius, max_radius, min_area, circularity_thresh, blur_size)
% 从二维数组中检测量子液滴
% 参数：
%     image_array: 二维数组（支持uint8或double类型，double类型会自动归一化）
%     density_thresh: 密度阈值（默认125）
%     min_radius: 最小检测半径（像素，默认1）
%     max_radius: 最大检测半径（像素，默认100）
%     min_area: 最小轮廓面积（像素，默认10）
%     circularity_thresh: 圆形度阈值（0-1，默认0.5）
%     blur_size: 高斯模糊核大小（奇数，默认5）
% 返回：
%     droplets: cell数组，每个元素为[x, y, radius]

% 设置默认参数
if nargin < 2 || isempty(density_thresh), density_thresh = 125; end
if nargin < 3 || isempty(min_radius), min_radius = 1; end
if nargin < 4 || isempty(max_radius), max_radius = 100; end
if nargin < 5 || isempty(min_area), min_area = 10; end
if nargin < 6 || isempty(circularity_thresh), circularity_thresh = 0.5; end
if nargin < 7 || isempty(blur_size), blur_size = 5; end

% 标准化图像数据
if ~isa(image_array, 'uint8')
    normalized = im2uint8(mat2gray(image_array));
else
    normalized = image_array;
end

% 预处理：高斯模糊
blurred = imgaussfilt(normalized, blur_size/3); % 注意：MATLAB的高斯模糊参数与OpenCV不同

% 二值化
thresh = blurred > density_thresh;

% 形态学操作去除噪声
kernel = strel('square', 3);
cleaned = imopen(thresh, kernel);
cleaned = imopen(cleaned, kernel); % 执行两次开运算

% 查找轮廓
[B, L] = bwboundaries(cleaned, 'noholes');

droplets = {};
for k = 1:length(B)
    boundary = B{k};
    
    % 计算面积
    area = polyarea(boundary(:,2), boundary(:,1));
    if area < min_area
        continue
    end
    
    % 计算周长和圆形度
    perimeter = arclength(boundary);
    if perimeter == 0
        continue
    end
    circularity = 4 * pi * area / (perimeter ^ 2);
    if circularity < circularity_thresh
        continue
    end
    
    % 最小包围圆
    [center, radius] = minEnclosingCircle(boundary);
    
    % 半径筛选
    if radius >= min_radius && radius <= max_radius
        droplets{end+1} = [round(center), round(radius)];
    end
end

    function perimeter = arclength(boundary)
        % 计算边界周长
        diff_coords = diff(boundary, 1);
        perimeter = sum(sqrt(sum(diff_coords.^2, 2)));
    end

    function [center, radius] = minEnclosingCircle(boundary)
        % 最小包围圆近似实现
        x = boundary(:,2);
        y = boundary(:,1);
        mx = mean(x);
        my = mean(y);
        dx = x - mx;
        dy = y - my;
        d = sqrt(dx.^2 + dy.^2);
        [max_d, idx] = max(d);
        center = [mx, my];
        radius = max_d;
        
        % 使用最远点作为初始圆心改进估计
        if max_d > 0
            center = [x(idx), y(idx)];
            d = sqrt((x - center(1)).^2 + (y - center(2)).^2);
            [max_d, ~] = max(d);
            radius = max_d;
        end
    end
end