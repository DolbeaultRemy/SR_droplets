function centers = calculateWeightedCenters(CC_split, imgCropped, selectFlag)
% 计算连通区域的加权中心
%
% 输入参数：
%   CC_split - bwconncomp函数返回的结构体，包含连通区域信息
%   imgCropped - 二维图像矩阵
%
% 输出参数：
%   centers - N×2矩阵，每行包含一个区域的加权中心坐标 [x, y]（列，行）

% 获取连通区域数量
numRegions = CC_split.NumObjects;

% 初始化中心坐标存储
centers = zeros(numRegions, 2);

% 遍历每个区域
for i = 1:numRegions

    % 获取当前区域的像素索引
    pixelIdx = CC_split.PixelIdxList{i};
    
    % 将线性索引转换为行列坐标
    [rows, cols] = ind2sub(size(imgCropped), pixelIdx);
    
    % 获取当前区域所有像素的值
    pixelValues = double(imgCropped(pixelIdx));  % 转换为double确保计算精度
    
    % 计算加权中心
    totalWeight = sum(pixelValues);
    
    if totalWeight > 0
        % 加权平均计算中心坐标
        centers(i, 1) = sum(cols .* pixelValues) / totalWeight;  % x坐标（列）
        centers(i, 2) = sum(rows .* pixelValues) / totalWeight;  % y坐标（行）
    else
        % 如果区域内所有像素值都为0，使用几何中心
        centers(i, 1) = mean(cols);
        centers(i, 2) = mean(rows);
        warning('区域 %d 的所有像素值为0，使用几何中心代替加权中心', i);
    end
end

centers = centers(selectFlag == 1, :);

end