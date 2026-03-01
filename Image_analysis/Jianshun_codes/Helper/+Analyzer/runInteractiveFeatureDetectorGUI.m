function runInteractiveFeatureDetectorGUI(od_imgs, scan_parameter_values, file_list, options, params)
%% runInteractiveFeatureDetectorGUI
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Interactive OD image viewer with patch detection
%   Supports slider, arrow keys, editable parameters, and displays scan parameter
%
% Notes:
%   Optional notes, references.

    numImages = numel(od_imgs);
    currentFrame = 1;

    %% --- Create figure ---
    % Try to find an existing figure by a unique tag
    hFig = findobj('Type','figure','Tag','InteractiveFeatureDetector');
    
    if isempty(hFig)
        % If not found, create a new figure
        hFig = figure('Name','OD Image Feature Detector', ...
                      'NumberTitle','off', ...
                      'Position',[100 100 1275 800], ...
                      'KeyPressFcn',@keyPressCallback, ...
                      'Tag','InteractiveFeatureDetector');  % <-- unique tag
    else
        % If figure exists, bring it to front
        figure(hFig);
        clf;
    end

    %% --- Axes ---
    hAx = axes('Parent',hFig,'Position',[0.025 0.1 0.6 0.85]);

    %% --- Frame slider ---
    hSlider = uicontrol('Style','slider','Min',1,'Max',numImages,'Value',1,...
        'SliderStep',[1/(numImages-1),10/(numImages-1)],...
        'Units','normalized','Position',[0.075 0.01 0.5 0.025],...
        'Callback',@(~,~) updateFrame());
    
    %% --- Parameter names ---
    paramNames = {'backgroundDiskFraction',...
                  'boundingBoxPadding',...
                  'dogGaussianSmallSigma',...
                  'dogGaussianLargeSigma',...
                  'adaptiveSensitivity', ...             %  .
                  'adaptiveNeighborhoodSize', ...              %   Larger → smoother masks, less sensitive to noise.
                  'minPeakFraction',...
                  'minimumPatchArea',...
                  'shapeMinArea', ...        
                  'shapeCloseRadius', ...
                  'shapeFillHoles', ...           
                  'intensityThreshFraction', ...
                  'edgeSigma', ...
                  'edgeThresholdLow', ...
                  'edgeThresholdHigh', ...
                  'pixelSize', ...
                  'magnification'};
    
    %% --- Parameter display names ---
    displayParamNames = {'Background Disk Fraction',...
                         'Bounding Box Padding (px)',...
                         'DoG Small Kernel Sigma',...
                         'DoG Large Kernel Sigma',...
                         'Adaptive threshold sensitivity',...
                         'Window size for adaptive threshold', ...
                         'Minimum Patch Peak Fraction',...
                         'Minimum Patch Area (px)',...
                         'Minimum Shape Area (px)', ...        
                         'Morphological Closing Radius (px)', ...
                         'Fill Internal Holes (0=false,1=true)', ...           
                         'Intensity Threshold Fraction', ...
                         'Canny Gaussian Smoothing Sigma', ...
                         'Canny Low Threshold Fraction', ...
                         'Canny High Threshold Fraction', ...
                         'Pixel Size (m/px)', ...
                         'Imaging System Magnification'};
    
    %% --- Parameter explanations ---
    paramDescriptions = {'Fraction of image used for background disk in morphological opening', ...
                         'Number of pixels to pad around detected cloud bounding box', ...
                         'Sigma of small Gaussian in Difference-of-Gaussians (DoG)', ...
                         'Sigma of large Gaussian in DoG', ...
                         'Global threshold that accounts for spatial variation of intensity', ...
                         'Local window size over which threshold is computed', ...
                         'Minimum fraction of max DoG response for detected patches', ...
                         'Minimum area (pixels) for detected patches', ...
                         'Minimum area (pixels) for internal shapes within patches', ...
                         'Radius (pixels) used for morphological closing to smooth shapes', ...
                         'Fill internal holes to ensure solid shapes (true/false)', ...
                         'Fraction of max intensity used to threshold features inside patches', ...
                         'Gaussian smoothing sigma used in Canny edge detection', ...
                         'Low threshold fraction for Canny edge detector', ...
                         'High threshold fraction for Canny edge detector', ...
                         'Physical size of one pixel (meters per pixel)', ...
                         'Magnification factor of imaging system'};

    
    nParams = numel(paramNames);
    hEdit = gobjects(nParams,1);
    
    for i = 1:nParams
        yTop = 0.9 - 0.05*(i-1);
        
        % Parameter name
        uicontrol('Style','text','Units','normalized',...
            'Position',[0.615 yTop 0.2 0.04],...
            'String',displayParamNames{i},'HorizontalAlignment','left', 'FontSize', 10, 'FontWeight','bold');
        
        % Parameter value edit box
        hEdit(i) = uicontrol('Style','edit','Units','normalized',...
            'Position',[0.875 yTop 0.1 0.04],...
            'String',num2str(params.(paramNames{i})),...
            'Callback',@(src,~) applyParams());
        
        % Explanation text below the edit box
        uicontrol('Style','text','Units','normalized',...
            'Position',[0.615 yTop - 0.01 0.35 0.03],...
            'String',paramDescriptions{i},...
            'HorizontalAlignment','left',...
            'FontSize',8,'ForegroundColor',[0.4 0.4 0.4], ...
            'BackgroundColor','none');
    end

    % Apply button
    uicontrol('Style','pushbutton','Units','normalized',...
        'Position',[0.7 yTop - 0.075 0.2 0.05],...
        'String','Apply',...
        'BackgroundColor',[0.2 0.6 1],...    % RGB values between 0 and 1
        'ForegroundColor','white',...         % Text color
        'FontSize',14, ...
        'FontWeight','bold',...
        'Callback',@(~,~) applyParams());


    % Initial plot
    updateFrame();

    %% --- Nested functions ---

    function applyParams()
        % Update params struct from textboxes
        for j = 1:nParams
            val = str2double(hEdit(j).String);
            if ~isnan(val)
                params.(paramNames{j}) = val;
            end
        end
        updateFrame();
    end

    function updateFrame()
        % Get current image
        idxImg   = round(get(hSlider,'Value'));
        img      = od_imgs{idxImg};  % numeric image
        
        [Ny, Nx] = size(img);

        dx       = params.pixelSize / params.magnification;
        dy       = dx;
        xAxis    = ((1:Nx)-(Nx+1)/2) * dx * 1e6; % μm
        yAxis    = ((1:Ny)-(Ny+1)/2) * dy * 1e6; % μm

        % Step 1: DoG detection
        [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart] = Analyzer.detectPatches(img, params);

        results                                                        = Analyzer.extractAndClassifyShapes(imgCropped, patchProps, xStart, yStart, params);
        
        %% Plotting in physical units (μm) ---
        cla(hAx);
        imagesc(hAx, xAxis, yAxis, img);
        axis(hAx,'equal','tight'); colormap(hAx, Colormaps.inferno());
        set(hAx,'FontSize',14,'YDir','normal');
        xlabel(hAx,'x (\mum)','FontSize',14); ylabel(hAx,'y (\mum)','FontSize',14);
        hold(hAx,'on');
        
        % Draw diagonal overlays
        Helper.drawODOverlays(xAxis(1), yAxis(1), xAxis(end), yAxis(end));
        Helper.drawODOverlays(xAxis(end), yAxis(1), xAxis(1), yAxis(end));
        
        % Plot patch centroids
        if ~isempty(patchCentroidsGlobal)
            patchCentroidsGlobal_um = [(patchCentroidsGlobal(:,1)-(Nx+1)/2)*dx*1e6, ...
                                       (patchCentroidsGlobal(:,2)-(Ny+1)/2)*dy*1e6];
            plot(hAx, patchCentroidsGlobal_um(:,1), patchCentroidsGlobal_um(:,2), 'ro','MarkerSize',4,'LineWidth',1);
        end
        
        % Plot patch ellipses
        for k = 1:numel(patchProps)
            a              = patchProps(k).MajorAxisLength/2 * dx * 1e6;
            b              = patchProps(k).MinorAxisLength/2 * dy * 1e6;
            phi            = deg2rad(-patchProps(k).Orientation);
            theta          = linspace(0,2*pi,100);
            R              = [cos(phi) -sin(phi); sin(phi) cos(phi)];
            ellipseXY      = [a*cos(theta(:)) b*sin(theta(:))]*R';
            cx_um          = ((patchProps(k).Centroid(1)+xStart-1)-(Nx+1)/2)*dx*1e6;
            cy_um          = ((patchProps(k).Centroid(2)+yStart-1)-(Ny+1)/2)*dy*1e6;
            ellipseXY(:,1) = ellipseXY(:,1) + cx_um;
            ellipseXY(:,2) = ellipseXY(:,2) + cy_um;
            plot(hAx, ellipseXY(:,1), ellipseXY(:,2),'g-','LineWidth',1);
        end
        
        % Overlay shapes
        for k = 1:numel(results)
            for n = 1:numel(results(k).boundaries)
                bnd = results(k).boundaries{n};
                bx  = ((bnd(:,2) - (Nx+1)/2)) * dx * 1e6;
                by  = ((bnd(:,1) - (Ny+1)/2)) * dy * 1e6;
                plot(hAx, bx, by, 'c-', 'LineWidth', 1.5);
            end
        end
        
        hold(hAx,'off');

        % Extract only filename (without path)
        [~, fname, ext]    = fileparts(file_list{idxImg});
        shortName          = [fname, ext];
        
        % Update figure title with shot + filename
        hFig.Name          = sprintf('Shot %d | %s', idxImg, shortName);

        % Update figure title with image index and patch count
        title(hAx, sprintf('Image %d/%d: Detected %d patches', ...
            idxImg, numImages, numel(patchProps)), ...
            'FontSize',16,'FontWeight','bold');
        
        % Text handle for scan parameter display
        txtHandle = text(hAx, 0.975, 0.975, '', ...
                    'Color', 'white', 'FontWeight', 'bold', ...
                    'FontSize', 24, 'Interpreter', 'tex', ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top');
        
        %% --- Generalized unit handling ---
        % Extract parameter row for this shot
        if iscell(scan_parameter_values)
            % Multi-parameter scan stored as cell array of row vectors
            param_row = scan_parameter_values{idxImg};
        else
            % Numeric vector / matrix
            param_row = scan_parameter_values(idxImg,:);
        end

        % Wrap single unit string into a cell if needed
        if ischar(options.scanParameterUnits) || isstring(options.scanParameterUnits)
            unitList = {char(options.scanParameterUnits)};
        else
            unitList = options.scanParameterUnits;  % assume cell array
        end

        % Ensure units list is long enough
        if numel(unitList) < numel(param_row)
            unitList(end+1:numel(param_row)) = {''};  % pad with empty units
        end

        % Build text lines for each parameter
        txtLines = cell(1, numel(param_row));
        for j = 1:numel(param_row)
            [unitSuffix, txtInterpreter] = getUnitInfo(unitList{j});
            txtLines{j} = sprintf('%.2f%s', param_row(j), unitSuffix);
        end

        % Join multiple parameters with newline
        txtHandle.String = strjoin(txtLines, '\n');
        txtHandle.Interpreter = txtInterpreter;  % use last parameter's interpreter
        
        % Update slider value
        hSlider.Value = idxImg;
        drawnow;
    end

    function keyPressCallback(event)
        switch event.Key
            case 'rightarrow'
                if currentFrame<numImages
                    currentFrame=currentFrame+1;
                    set(hSlider,'Value',currentFrame);
                    updateFrame();
                end
            case 'leftarrow'
                if currentFrame>1
                    currentFrame=currentFrame-1;
                    set(hSlider,'Value',currentFrame);
                    updateFrame();
                end
        end
    end
end

%% --- Helper function ---
function [unitSuffix, txtInterpreter] = getUnitInfo(u)
    switch lower(u)
        case {'degrees','deg','°'}
            unitSuffix     = '^\circ';
            txtInterpreter = 'tex';
        case {'gauss','g'}
            unitSuffix     = ' G';
            txtInterpreter = 'none';
        otherwise
            unitSuffix     = '';
            txtInterpreter = 'none';
    end
end