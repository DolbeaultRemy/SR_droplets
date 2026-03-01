function plotMeanWithSE(scan_values, data_values, varargin)
%% plotMeanWithSE
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots mean ± standard error vs a scan parameter.
%
% Notes:
%   Optional notes, references.

    % --- Parse optional name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'YLim', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'mean_with_se.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'HoldOn', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Compute mean and standard error ---
    [unique_vals, ~, idx] = unique(scan_values);
    mean_vals = zeros(size(unique_vals));
    stderr_vals = zeros(size(unique_vals));
    for i = 1:length(unique_vals)
        if iscell(data_values)
            group = data_values{idx == i};
        else
            group = data_values(idx == i);
        end
        if iscell(group)
            groupVals = [group{:}];
        else
            groupVals = group;
        end
        mean_vals(i)   = mean(groupVals);
        stderr_vals(i) = std(groupVals) / sqrt(length(groupVals));
    end

    % --- Create figure ---
    if isempty(opts.FigNum)
        fig = figure;
    else
        fig = figure(opts.FigNum);
    end
    if ~opts.HoldOn
        clf(fig);
    end
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);

    % --- Plot error bars ---
    errorbar(unique_vals, mean_vals, stderr_vals, 'o--', ...
        'LineWidth', 1.8, 'MarkerSize', 6, 'CapSize', 5);

    if opts.HoldOn
        hold on;
    end

    % --- Axis formatting ---
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~isempty(opts.YLim)
        ylim(opts.YLim);
    end
    xlabel(opts.XLabel, 'Interpreter', 'latex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel(opts.YLabel, 'Interpreter', 'latex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title(opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2, 'FontWeight', 'bold');
    grid on;

    % --- Save figure ---
    Plotter.saveFigure(fig, ...
        'SaveFileName', opts.SaveFileName, ...
        'SaveDirectory', opts.SaveDirectory, ...
        'SkipSaveFigures', opts.SkipSaveFigures);

end