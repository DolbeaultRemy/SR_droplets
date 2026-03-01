function runInteractiveODImageViewer(od_imgs, scan_parameter_values, file_list, options)
%% runInteractiveFeatureDetectorGUI
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Interactive OD image viewer with patch detection
%   Supports slider, arrow keys, editable parameters, and displays scan parameter
%
% Inputs:
%   od_imgs                : cell array of 2D OD images
%   scan_parameter_values  : array of corresponding scan parameter values
%   file_list              : cell array of corresponding filenames
%   options                : struct with fields
%       .pixel_size, .magnification, .center, .font, .zoom_size, .scan_parameter
%
% Notes:
%   Optional notes, references.

    %% --- Create Figure ---
    % Try to find an existing figure by a unique tag
    hFig = findobj('Type','figure','Tag','InteractiveImageViewer');
    
    if isempty(hFig)
        % If not found, create a new figure
        hFig = figure('Name','OD Image Viewer', ...
                      'NumberTitle','off', ...
                      'Position', [50 50 1000 800], ...
                      'KeyPressFcn',@keyPressCallback, ...
                      'Tag','InteractiveImageViewer');  % <-- unique tag
    else
        % If figure exists, bring it to front
        figure(hFig);
        clf;
    end

    %% --- Get image size ---
    [Ny, Nx] = size(od_imgs{1});
    
    %% --- Pixel size and axes in μm ---
    dx = options.pixel_size / options.magnification;
    dy = dx;  % square pixels
    x = ((1:Nx) - (Nx+1)/2) * dx * 1e6;  
    y = ((1:Ny) - (Ny+1)/2) * dy * 1e6;
    
    %% --- Display first image ---
    hAx = axes('Parent', hFig);
    hImg = imagesc(hAx, x, y, od_imgs{1});
    axis(hAx, 'equal', 'tight')
    colormap(hAx, Colormaps.inferno());
    set(hAx, 'FontSize', 14, 'YDir', 'normal');
    xlabel(hAx, 'x (\mum)', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', options.font);
    ylabel(hAx, 'y (\mum)', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', options.font);
    title(hAx, ['Measurement: ', options.titleString], 'FontSize', 16, ...
     'FontWeight', 'bold', 'Interpreter', 'tex', 'FontName', options.font);
    colorbarHandle = colorbar(hAx);
    ylabel(colorbarHandle, 'Optical Density', 'Rotation', -90, 'FontSize', 14, 'FontName', options.font);

    hold(hAx, 'on')
    % Draw diagonal overlays once
    Helper.drawODOverlays(x(1), y(1), x(end), y(end));
    Helper.drawODOverlays(x(end), y(1), x(1), y(end));
    hold(hAx, 'off')
    
    txtHandle = text(hAx, 0.975, 0.975, '', ...
        'Color', 'white', 'FontWeight', 'bold', ...
        'FontSize', 24, 'Interpreter', 'tex', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top');

    % Slider
    sliderHandle = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', length(od_imgs), 'Value', 1, ...
        'SliderStep', [1/(length(od_imgs)-1), 10/(length(od_imgs)-1)], ...
        'Position', [150 5 700 20], ...
        'Callback', @(src, ~) updateImage(round(src.Value)));

    % Initialize
    currentIdx = 1;
    updateImage(currentIdx);

    % Arrow key callback
    set(hFig, 'KeyPressFcn', @(src, event) keyPressCallback(event));

    %% --- Nested Functions ---
    function updateImage(idx)
        currentIdx = idx;
        hImg.CData = od_imgs{idx};

        % Extract only filename (without path)
        [~, fname, ext] = fileparts(file_list{idx});
        shortName = [fname, ext];

        % Update figure title with shot + filename
        hFig.Name = sprintf('Shot %d | %s', idx, shortName);

        %% --- Generalized unit handling ---
        % Extract parameter row for this shot
        if iscell(scan_parameter_values)
            % Multi-parameter scan stored as cell array of row vectors
            param_row = scan_parameter_values{idx};
        else
            % Numeric vector / matrix
            param_row = scan_parameter_values(idx,:);
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

        sliderHandle.Value = idx;
        drawnow;
    end

    function keyPressCallback(event)
        switch event.Key
            case 'rightarrow'
                if currentIdx < length(od_imgs)
                    updateImage(currentIdx + 1);
                end
            case 'leftarrow'
                if currentIdx > 1
                    updateImage(currentIdx - 1);
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