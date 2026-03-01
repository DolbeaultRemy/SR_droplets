function plotHeatmap(results, scan_parameter_values, scan_reference_values, fieldName, varargin)
%% plotHeatmap
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots a heatmap for a field in a struct array.
%
% Notes:
%   Optional notes, references.

    % --- Parse optional inputs ---
    p = inputParser;
    addParameter(p, 'FigNum', []);
    addParameter(p, 'Colormap', parula);
    addParameter(p, 'CLim', []);
    addParameter(p, 'ColorScale', 'linear', @(x) any(validatestring(x,{'linear','log'}))); 
    addParameter(p, 'XLabel', 'Param1');
    addParameter(p, 'YLabel', 'Param2');
    addParameter(p, 'Title', fieldName);
    addParameter(p, 'FontName', 'Arial');
    addParameter(p, 'FontSize', 14);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'heatmap.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'PlotDirectly', false, @islogical);  
    parse(p, varargin{:});
    opts = p.Results;

    % --- Convert reference/parameter values to numeric arrays ---
    scanRefArray   = cell2mat(scan_reference_values(:));   % nParams × 2
    scanParamArray = cell2mat(scan_parameter_values(:));   % nTotal × 2

    % --- Determine unique X and Y values ---
    XValues = sort(unique(scanRefArray(:,2)), 'ascend'); % Param1 → X
    YValues = sort(unique(scanRefArray(:,1)), 'ascend'); % Param2 → Y
    nX = numel(XValues);
    nY = numel(YValues);

    % --- Preallocate data matrix ---
    data_matrix = NaN(nY, nX);

    % --- Fill matrix by looping over unique parameter pairs ---
    for k = 1:size(scanRefArray,1)
        X = scanRefArray(k,2);
        Y = scanRefArray(k,1);

        row = find(ismembertol(YValues, Y, 1e-6));
        col = find(ismembertol(XValues, X, 1e-6));
        if isempty(row) || isempty(col), continue; end

        repIdx = find(ismembertol(scanParamArray, [Y X], 1e-6, 'ByRows', true));
        if isempty(repIdx), continue; end

        if opts.PlotDirectly
            data_matrix(row,col) = results.(fieldName)(k);
        else
            vals = results.(fieldName)(repIdx);
            data_matrix(row,col) = mean(vals, 'omitnan');
        end
    end

    % --- Create figure ---
    if isempty(opts.FigNum)
        fig = figure;
    else
        fig = figure(opts.FigNum);
    end
    clf(fig);
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);

    % --- Plot heatmap ---
    imagesc(XValues, YValues, data_matrix);
    colormap(feval(opts.Colormap));
    if ~isempty(opts.CLim)
        caxis(opts.CLim);
    end
    set(gca, ...
        'FontName', opts.FontName, ...
        'FontSize', opts.FontSize, ...
        'YDir', 'normal', ...
        'ColorScale', opts.ColorScale);

    % --- Labels and title ---
    xlabel(opts.XLabel, 'Interpreter', 'tex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel(opts.YLabel, 'Interpreter', 'tex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title(opts.Title, 'Interpreter', 'latex', 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2, 'FontWeight', 'bold');
    colorbar;

    % --- Save figure ---
    Plotter.saveFigure(fig, ...
        'SaveFileName', opts.SaveFileName, ...
        'SaveDirectory', opts.SaveDirectory, ...
        'SkipSaveFigures', opts.SkipSaveFigures);
end
