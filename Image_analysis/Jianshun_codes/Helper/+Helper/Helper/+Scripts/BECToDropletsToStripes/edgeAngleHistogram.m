% edgeAngleHistogram.m
function [allAngles, edgeLengths, edgesPlot] = edgeAngleHistogram(points)
    % 计算边之间夹角的直方图
    % 输入: points - 点坐标矩阵 (N×2)
    % 输出: allAngles - 所有计算得到的夹角（度数）
    %       edgeLengths - 边的长度
    %       edgesPlot - 边的索引
    
    % 检查输入有效性
    if isempty(points) || size(points, 1) < 3
        allAngles = [];
        edgeLengths = [];
        edgesPlot = [];
        fprintf('点数不足，无法进行三角剖分\n');
        return;
    end
    
    % Delaunay三角剖分
    DT = delaunayTriangulation(points);
    
    % 获取凸包上的点（外层点）
    convexHullIdx = convexHull(DT);
    isOuterPoint = false(size(points, 1), 1);
    isOuterPoint(convexHullIdx) = true;
    
    % 获取所有边
    edges = DT.edges;
    edgesPlot = [];
    edgeLengths = [];
    edgeInner = [];
    
    % 选择边（排除连接两个外层点的边）
    for i = 1:size(edges, 1)
        edge = edges(i, :);
        p1 = points(edge(1), :);
        p2 = points(edge(2), :);
        
        % 如果这条边不是连接两个外层点的，则保留
        if ~(isOuterPoint(edge(1)) && isOuterPoint(edge(2)))
            edgeLengths = [edgeLengths; norm(p1 - p2)];
            edgeInner = [edgeInner; edge(1), edge(2)];
        end
    end
    
    % 如果没有有效的边，返回空
    if isempty(edgeInner)
        allAngles = [];
        return;
    end
    
    % 计算每条边的向量表示
    vectors = zeros(size(edgeInner, 1), 2);
    for i = 1:size(edgeInner, 1)
        p1 = points(edgeInner(i, 1), :);
        p2 = points(edgeInner(i, 2), :);
        vectors(i, :) = p2 - p1;  % 边向量
    end
    
    % 计算所有边对之间的夹角
    allAngles = calculate_angles_between_edges(vectors);
    edgesPlot = edgeInner;
end

function allAngles = calculate_angles_between_edges(vectors)
    % 计算边向量之间的夹角
    
    numEdges = size(vectors, 1);
    allAngles = [];
    
    if numEdges < 2
        fprintf('边数不足，无法计算夹角\n');
        return;
    end
    
    % 计算每对边之间的夹角
    for i = 1:numEdges-1
        for j = i+1:numEdges
            % 获取两个向量
            v1 = vectors(i, :);
            v2 = vectors(j, :);
            
            % 归一化向量
            norm_v1 = norm(v1);
            norm_v2 = norm(v2);
            
            if norm_v1 > 0 && norm_v2 > 0
                % 计算点积
                dot_product = dot(v1/norm_v1, v2/norm_v2);
                
                % 限制点积范围，防止数值误差
                dot_product = max(min(dot_product, 1), -1);
                
                % 计算夹角（弧度）并转换为度数
                angle_rad = acos(dot_product);
                angle_deg = angle_rad * 180 / pi;
                
                % 添加到夹角列表
                allAngles = [allAngles; angle_deg];
            end
        end
    end
end