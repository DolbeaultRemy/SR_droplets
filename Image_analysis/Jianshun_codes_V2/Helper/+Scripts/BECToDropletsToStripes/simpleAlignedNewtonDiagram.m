function [allRotatedPoints, edgeLengths, edgesPlot] = simpleAlignedNewtonDiagram(points, N)
    % 计算对齐的牛顿图
    % 输入:
    %   points - 点坐标矩阵 (N×2)
    %   N - (可选) 选择第N短的边进行对齐，默认值为1（最短的边）
    %       当总边数小于N时，使用最长的边
    % 输出:
    %   allRotatedPoints - 旋转后的点坐标
    %   edgeLengths - 边的长度
    %   edgesPlot - 边的索引
    
    % 设置默认参数
    if nargin < 2
        N = 3;  % 默认为最短边
    end
    
    % 检查输入有效性
    if isempty(points) || size(points, 1) < 3
        allRotatedPoints = [];
        edgeLengths = [];
        edgesPlot = [];
        return;
    end
    
    % Delaunay三角剖分
    DT = delaunayTriangulation(points);
    
    % 获取凸包上的点（外层点）
    convexHullIdx = convexHull(DT);
    isOuterPoint = false(size(points, 1), 1);
    isOuterPoint(convexHullIdx) = true;
    innerPoints = points(~isOuterPoint, :);
    innerPointIndices = find(~isOuterPoint);
    
    % 获取所有边
    edges = DT.edges;
    edgesPlot = [];
    
    % 选择边（排除连接两个外层点的边）
    edgeLengths = [];
    edgeInner = [];
    
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
    
    % 如果没有内层点，处理特殊情况
    if isempty(innerPoints)
        [edgeLengths, edgeInner] = handle_no_inner_points(points, edges);
    end
    
    % 处理每个内层点，计算旋转后的邻居点
    allRotatedPoints = [];
    numInner = size(innerPoints, 1);
    
    for i = 1:numInner
        centerPointIdx = innerPointIndices(i);
        centerPoint = points(centerPointIdx, :);
        
        % 找到与当前点相连的边
        connectedEdges = edgeInner(any(edgeInner == centerPointIdx, 2), :);
        
        % 获取所有邻居点索引（排除自身）
        neighborIndices = unique(connectedEdges(connectedEdges ~= centerPointIdx));
        
        if isempty(neighborIndices)
            continue;
        end
        
        % 计算旋转后的邻居点（使用第N短的边作为参考）
        rotatedNeighbors = rotate_neighbors_to_reference(centerPoint, points, neighborIndices, N);
        allRotatedPoints = [allRotatedPoints; rotatedNeighbors];
    end
    
    edgesPlot = edgeInner;
end

function rotatedNeighbors = rotate_neighbors_to_reference(centerPoint, points, neighborIndices, N)
    % 将邻居点旋转到以第N短的边为参考的方向
    % 如果总邻居数小于N，则使用最长的边
    
    % 获取邻居点
    neighborPoints = points(neighborIndices, :);
    
    % 计算到每个邻居的距离
    distances = sqrt(sum((neighborPoints - centerPoint).^2, 2));
    
    % 对距离进行排序
    [sortedDistances, sortedIndices] = sort(distances);
    
    % 选择参考边
    if N <= length(sortedDistances)
        % 选择第N短的边
        refIdx = sortedIndices(N);
    else
        % 如果总邻居数小于N，选择最长的边
        refIdx = sortedIndices(end);
    end
    
    % 获取参考点
    referencePoint = neighborPoints(refIdx, :);
    
    % 移除参考点，保留其他邻居点
    neighborPoints(refIdx, :) = [];
    
    % 如果没有其他邻居点，返回空
    if isempty(neighborPoints)
        rotatedNeighbors = [];
        return;
    end
    
    % 计算旋转角度（使参考点位于0度方向）
    vectorToReference = referencePoint - centerPoint;
    rotationAngle = -atan2(vectorToReference(2), vectorToReference(1));
    
    % 旋转所有邻居点
    rotationMatrix = [cos(rotationAngle), -sin(rotationAngle);
                     sin(rotationAngle), cos(rotationAngle)];
    
    % 计算邻居点相对于中心点的坐标并旋转
    neighborVectors = neighborPoints - centerPoint;
    rotatedNeighbors = (rotationMatrix * neighborVectors')';
end

function [edgeLengths, edgeInner] = handle_no_inner_points(points, edges)
    % 处理没有内层点的情况：为每个点保留最短的边
    
    numPoints = size(points, 1);
    shortestEdgeInfo = cell(numPoints, 1);
    edgeLengths = [];
    edgeInner = [];
    
    % 第一遍：找到每个点的最短边
    for i = 1:size(edges, 1)
        edge = edges(i, :);
        point1 = edge(1);
        point2 = edge(2);
        p1 = points(point1, :);
        p2 = points(point2, :);
        currentLength = norm(p1 - p2);
        
        % 更新点1的最短边
        if isempty(shortestEdgeInfo{point1}) || currentLength < shortestEdgeInfo{point1}.length
            shortestEdgeInfo{point1} = struct('edge', edge, 'length', currentLength);
        end
        
        % 更新点2的最短边
        if isempty(shortestEdgeInfo{point2}) || currentLength < shortestEdgeInfo{point2}.length
            shortestEdgeInfo{point2} = struct('edge', edge, 'length', currentLength);
        end
    end
    
    % 第二遍：收集所有最短边（去重）
    collectedEdges = [];
    for i = 1:numPoints
        if ~isempty(shortestEdgeInfo{i})
            currentEdge = sort(shortestEdgeInfo{i}.edge);
            
            % 检查是否已经收集过这条边
            isNew = true;
            for j = 1:size(collectedEdges, 1)
                if all(collectedEdges(j, :) == currentEdge)
                    isNew = false;
                    break;
                end
            end
            
            if isNew
                collectedEdges = [collectedEdges; currentEdge];
                edgeLengths = [edgeLengths; shortestEdgeInfo{i}.length];
                edgeInner = [edgeInner; shortestEdgeInfo{i}.edge];
            end
        end
    end
end