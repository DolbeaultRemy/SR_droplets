function plotG2PDF(results, theta_query, varargin)
%% plotG2PDF
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots 2D heatmap of PDFs of g²(θ) values for different scan parameters.
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
    addParameter(p, 'SaveFileName', 'G2PDF.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'NumPoints', 200, @(x) isscalar(x));
    addParameter(p, 'DataRange', [], @(x) isempty(x) || numel(x)==2);
    addParameter(p, 'XLim', [], @(x) isempty(x) || numel(x)==2);
    addParameter(p, 'Colormap', @jet);
    addParameter(p, 'ColorScale', 'linear', @(x) any(validatestring(x,{'linear','log'})));
    addParameter(p, 'PlotType', 'histogram', @(x) any(validatestring(x,{'kde','histogram'})));
    addParameter(p, 'NumberOfBins', 50, @isscalar);
    addParameter(p, 'NormalizeHistogram', true, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Extract data at requested theta ---
    theta_vals = results.theta_values;
    [~, idx_theta] = min(abs(theta_vals - theta_query)); % closest index
    N_params = numel(results.g2_curves);
    dataCell = cell(N_params,1);
    for i = 1:N_params
        G = results.g2_curves{i};      % [N_reps × Nθ]
        dataCell{i} = G(:, idx_theta);  % all repetitions at chosen theta
    end

    referenceValues = results.scan_parameter_values;

    % --- Determine y-range ---
    if isempty(opts.DataRange)
        allData = cell2mat(dataCell(:));
        y_min   = min(allData);
        y_max   = max(allData);
    else
        y_min   = opts.DataRange(1);
        y_max   = opts.DataRange(2);
    end

    % --- Prepare PDF grid/matrix ---
    if strcmpi(opts.PlotType,'kde')
        y_grid     = linspace(y_min, y_max, opts.NumPoints);
        pdf_matrix = zeros(numel(y_grid), N_params);
    else
        edges      = linspace(y_min, y_max, opts.NumberOfBins+1);
        binCenters = (edges(1:end-1) + edges(2:end))/2;
        pdf_matrix = zeros(numel(binCenters), N_params);
    end

    % --- Compute PDFs ---
    for i = 1:N_params
        data = dataCell{i};
        data = data(~isnan(data));
        if isempty(data), continue; end

        if strcmpi(opts.PlotType,'kde')
            f = ksdensity(data, y_grid);
            pdf_matrix(:,i) = f;
        else
            counts = histcounts(data, edges);
            if opts.NormalizeHistogram
                binWidth = edges(2) - edges(1);
                counts   = counts / (sum(counts) * binWidth);
            end
            pdf_matrix(:,i) = counts(:);
        end
    end

    % --- Plot heatmap ---
    fig = figure(opts.FigNum); clf(fig);
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);
    if strcmpi(opts.PlotType,'kde')
        imagesc(referenceValues, y_grid, pdf_matrix);
    else
        imagesc(referenceValues, binCenters, pdf_matrix);
    end

    set(gca, 'YDir', 'normal', 'FontName', opts.FontName, 'FontSize', opts.FontSize, 'ColorScale', opts.ColorScale);
    xlabel(opts.XLabel, 'Interpreter', 'tex', 'FontSize', opts.FontSize, 'FontName', opts.FontName);
    ylabel(opts.YLabel, 'Interpreter', 'tex', 'FontSize', opts.FontSize, 'FontName', opts.FontName);

    % --- Compute fraction of pi ---
    theta_factor = theta_query / pi;  % e.g., pi/6 -> 1/6
    
    % --- Format as string ---
    if abs(round(1/theta_factor) - 1/theta_factor) < 1e-6  % e.g., 1/6
        theta_str = sprintf('\\pi/%d', round(1/theta_factor));
    elseif abs(theta_factor - round(theta_factor)) < 1e-6   % e.g., pi or 2*pi
        theta_str = sprintf('%d\\pi', round(theta_factor));
    else
        theta_str = sprintf('%.3g \\pi', theta_factor);  % fallback numeric
    end
    
    % --- Set title ---
    title(sprintf('%s | g^{(2)}(\\delta\\theta = %s)', opts.Title, theta_str), ...
          'FontName', opts.FontName, 'FontSize', opts.FontSize + 2, 'FontWeight', 'bold');
    
    colormap(feval(opts.Colormap));
    c = colorbar;
    if strcmpi(opts.PlotType,'kde')
        ylabel(c, 'PDF', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    else
        if opts.NormalizeHistogram
            ylabel(c, 'Probability Density', 'Rotation', -90, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        else
            ylabel(c, 'Counts', 'Rotation', -90, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        end
    end

    if ~isempty(opts.XLim)
        xlim(opts.XLim);
    end

    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end
