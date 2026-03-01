% calculate_aspect_ratios_simple.m
function [aspect_ratio] = calculate_aspect_ratios_simple(CC_split)
    % 计算分割区域的长宽比
    % 输入: CC_split - 区域分割的结构体
    % 输出: aspect_ratio - 每个区域的长宽比数组
    
    % 获取区域属性，包括方向
    stats = regionprops(CC_split, 'BoundingBox', 'Centroid', 'Area', ...
                        'Orientation', 'MajorAxisLength', 'MinorAxisLength');
    
    num_regions = length(stats);
    aspect_ratio = zeros(num_regions, 1);
    
    % 计算每个区域的长宽比
    for i = 1:num_regions
        major_length = stats(i).MajorAxisLength;
        minor_length = stats(i).MinorAxisLength;
        aspect_ratio(i) = major_length / minor_length;
    end
end