function plotPDF(dataCell, referenceValues, varargin)
%% plotPDF
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots 2D heatmap of PDFs for grouped data (Histogram or KDE).
%
% Notes:
%   Optional notes, references.

    % --- Parse optional inputs ---
    p = inputParser;
    addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FigNum', 1, @(x) isscalar(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'pdf.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'NumPoints', 200, @(x) isscalar(x));
    addParameter(p, 'DataRange', [], @(x) isempty(x) || numel(x)==2);
    addParameter(p, 'XLim', [], @(x) isempty(x) || numel(x)==2);
    addParameter(p, 'Colormap', @jet);
    addParameter(p, 'PlotType', 'histogram', @(x) any(validatestring(x,{'kde','histogram'})));
    addParameter(p, 'NumberOfBins', 50, @isscalar);
    addParameter(p, 'NormalizeHistogram', true, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    N_params = numel(referenceValues);

    % --- Determine y-range ---
    if isempty(opts.DataRange)
        allData = cell2mat(dataCell(:));
        y_min   = min(allData);
        y_max   = max(allData);
    else
        y_min   = opts.DataRange(1);
        y_max   = opts.DataRange(2);
    end

    if strcmpi(opts.PlotType,'kde') % KDE
        y_grid     = linspace(y_min, y_max, opts.NumPoints);
        pdf_matrix = zeros(numel(y_grid), N_params);
    else % Histogram
        edges      = linspace(y_min, y_max, opts.NumberOfBins+1);
        binCenters = (edges(1:end-1) + edges(2:end)) / 2;
        pdf_matrix = zeros(numel(binCenters), N_params);
    end

    % --- Compute PDFs ---
    for i = 1:N_params
        data = dataCell{i};
        data = data(~isnan(data));
        if isempty(data), continue; end

        if strcmpi(opts.PlotType,'kde') % KDE
            f                = ksdensity(data, y_grid);
            pdf_matrix(:, i) = f;
        else % Histogram
            counts = histcounts(data, edges);
            if opts.NormalizeHistogram
                binWidth = edges(2) - edges(1);
                counts   = counts / (sum(counts) * binWidth); % probability density
            end
            pdf_matrix(:, i) = counts(:);
        end
    end

    % --- Plot heatmap ---
    fig = figure(opts.FigNum); clf(fig);
    set(fig, 'Color', 'w', 'Position',[100 100 950 750]);

    if strcmpi(opts.PlotType,'kde')
        imagesc(referenceValues, y_grid, pdf_matrix);
    else
        imagesc(referenceValues, binCenters, pdf_matrix);
    end

    set(gca, 'YDir', 'normal', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    xlabel(opts.XLabel, 'Interpreter', 'latex', 'FontSize', opts.FontSize, 'FontName', opts.FontName);
    ylabel(opts.YLabel, 'Interpreter', 'latex', 'FontSize', opts.FontSize, 'FontName', opts.FontName);
    title(opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2, 'FontWeight', 'bold');
    colormap(feval(opts.Colormap));
    c = colorbar;
    if strcmpi(opts.PlotType,'kde')
        ylabel(c, 'PDF', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    else
        if opts.NormalizeHistogram
            ylabel(c, 'Probability Density', 'Rotation', -90, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        else
            ylabel(c, 'Counts', 'FontName', 'Rotation', -90, opts.FontName, 'FontSize', opts.FontSize);
        end
    end

    if ~isempty(opts.XLim)
        xlim(opts.XLim);
    end

    % --- Save figure ---
    Plotter.saveFigure(fig, ...
        'SaveFileName', opts.SaveFileName, ...
        'SaveDirectory', opts.SaveDirectory, ...
        'SkipSaveFigures', opts.SkipSaveFigures);

end