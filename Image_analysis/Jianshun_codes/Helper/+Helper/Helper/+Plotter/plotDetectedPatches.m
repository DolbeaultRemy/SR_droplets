function plotDetectedPatches(img, patchProps, results, params, xStart, yStart, figTitle)
%% plotDetectedPatches
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots the image, DoG-detected patch ellipses, and extracted shape boundaries
%   with optional text overlay for Dice score and patch count.
%
% Inputs:
%   img        - original 2D image
%   patchProps - struct array from detectPatches
%   results    - struct array from extractAndClassifyShapes
%   params     - parameter struct
%   xStart, yStart - offsets for cropped regions
%   figTitle   - title string (e.g., 'BayesOpt Candidate')
%
% Notes:
%   Optional notes, references. 

    if nargin < 7
        figTitle = '';
    end

    [Ny, Nx] = size(img);
    dx = params.pixelSize / params.magnification;
    dy = dx;

    % Physical axes (μm)
    xAxis = ((1:Nx)-(Nx+1)/2) * dx * 1e6;
    yAxis = ((1:Ny)-(Ny+1)/2) * dy * 1e6;

    %% --- Create Figure ---
    % Try to find an existing figure by a unique tag
    hFig = findobj('Type','figure','Tag','OptimizerViewer');
    
    if isempty(hFig)
        % If not found, create a new figure
        hFig = figure('Name','Optimization Live Viewer', ...
                      'NumberTitle','off', ...
                      'Position', [200 200 600 600], ...
                      'KeyPressFcn',@keyPressCallback, ...
                      'Tag','OptimizerViewer');  % <-- unique tag
    else
        % If figure exists, bring it to front
        figure(hFig);
        clf;
    end

    imagesc(xAxis, yAxis, img); 
    axis(gca,'equal','tight'); 
    colormap(Colormaps.inferno());
    set(gca,'FontSize',14,'YDir','normal');
    xlabel(gca,'x (\mum)','FontName', 'Bahnschrift', 'FontSize', 14, 'FontWeight', 'bold'); 
    ylabel(gca,'y (\mum)','FontName', 'Bahnschrift', 'FontSize', 14, 'FontWeight', 'bold');
    
    hold on;

    % Draw patch ellipses
    for k = 1:numel(patchProps)
        a = patchProps(k).MajorAxisLength/2 * dx * 1e6;
        b = patchProps(k).MinorAxisLength/2 * dy * 1e6;
        phi = deg2rad(-patchProps(k).Orientation);
        theta = linspace(0,2*pi,100);
        R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
        ellipseXY = [a*cos(theta(:)) b*sin(theta(:))]*R';
        cx_um = ((patchProps(k).Centroid(1)+xStart-1)-(Nx+1)/2)*dx*1e6;
        cy_um = ((patchProps(k).Centroid(2)+yStart-1)-(Ny+1)/2)*dy*1e6;
        ellipseXY(:,1) = ellipseXY(:,1) + cx_um;
        ellipseXY(:,2) = ellipseXY(:,2) + cy_um;
        plot(ellipseXY(:,1), ellipseXY(:,2),'g-','LineWidth',1);
    end

    % Overlay shapes
    for k = 1:numel(results)
        for n = 1:numel(results(k).boundaries)
            bnd = results(k).boundaries{n};
            bx = ((bnd(:,2)-(Nx+1)/2)) * dx * 1e6;
            by = ((bnd(:,1)-(Ny+1)/2)) * dy * 1e6;
            plot(bx, by, 'c-', 'LineWidth', 1.5);
        end
    end

    % Compute Dice if gtMask is available
    if isfield(params,'gtMask') && ~isempty(params.gtMask)
        predictedMask = false(size(img));
        for k = 1:numel(results)
            BW = results(k).BW;
            if isempty(BW), continue; end
            [h,w] = size(BW);
            yIdx = (yStart:yStart+h-1);
            xIdx = (xStart:xStart+w-1);
            predictedMask(yIdx, xIdx) = predictedMask(yIdx, xIdx) | BW;
        end
        diceScore = dice(predictedMask, params.gtMask);
    
        % Convert logical mask to coordinates for overlay
        [maskY, maskX] = find(params.gtMask);
        maskX_um = (maskX - (Nx+1)/2) * dx * 1e6;
        maskY_um = (maskY - (Ny+1)/2) * dy * 1e6;
    
        % GT mask Overlay 
        overlayColor = [0.94 0.94 0.94];  % RGB from middle of coolwarm
        scatter(maskX_um, maskY_um, 15, overlayColor, 'filled', ...
                'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    else
        diceScore = NaN;
    end

    % Add text overlay in top-right
    txtStr = sprintf('Detected patches: %d\nDice score: %.3f', numel(patchProps), diceScore);
    text(0.975, 0.975, txtStr, 'Units','normalized', 'Color','w', ...
        'FontName', 'Bahnschrift', 'FontWeight','bold', 'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','top');

    title(figTitle,'FontName', 'Bahnschrift', 'FontSize',16,'FontWeight','bold');
    hold off;
end
