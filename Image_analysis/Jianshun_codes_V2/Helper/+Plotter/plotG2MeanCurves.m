function plotG2MeanCurves(g2_mean_all, g2_error_all, theta_values, scan_parameter_values, scan_parameter_units, varargin)
%% plotG2MeanCurves
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots mean g2 angular correlations with optional parameters.
%
% Notes:
%   Optional notes, references.

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', 'g^{(2)}(\delta\theta) vs \delta\theta', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', '$\delta\theta / \pi$', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLabel', '$g^{(2)}(\delta\theta)$', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'Colormap', @parula);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'figure.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'YLim', [0 1], @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    nParams = size(g2_mean_all, 2);
    
    % Determine unit suffix and interpreter
    switch lower(scan_parameter_units)
        case {'degrees', 'deg', '°'}
            unitSuffix = '^\circ';   % LaTeX degree symbol
            txtInterpreter = 'tex';
        case {'gauss', 'g'}
            unitSuffix = ' G';
            txtInterpreter = 'none';
        otherwise
            unitSuffix = '';
            txtInterpreter = 'none';
    end

    % --- Create figure ---
    if isempty(opts.FigNum)
        fig = figure;
    else
        fig = figure(opts.FigNum);
    end
    clf(fig);
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);
    hold on;

    % --- Colormap ---
    cmap = opts.Colormap(nParams);
    
    % --- Plot data with errorbars ---
    legend_entries = cell(nParams, 1);
    for i = 1:nParams
        errorbar(theta_values/pi, g2_mean_all{i}, g2_error_all{i}, ...
            'o', 'Color', cmap(i,:), 'MarkerSize', 4, 'MarkerFaceColor', cmap(i,:), 'CapSize', 4);

        % Update overlay text with scan parameter + unit
        legend_entries{i} = sprintf('%.2f%s', scan_parameter_values(i), unitSuffix);
    end

    % --- Formatting ---
    xlabel(opts.XLabel, 'Interpreter', 'latex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel(opts.YLabel, 'Interpreter', 'latex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title(opts.Title, 'FontName', opts.FontName, 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2);
    legend(legend_entries, 'Interpreter', txtInterpreter, 'Location', 'bestoutside');
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylim(opts.YLim);
    grid on;

    % --- Save figure ---
    Plotter.saveFigure(fig, ...
        'SaveFileName', opts.SaveFileName, ...
        'SaveDirectory', opts.SaveDirectory, ...
        'SkipSaveFigures', opts.SkipSaveFigures);

end
