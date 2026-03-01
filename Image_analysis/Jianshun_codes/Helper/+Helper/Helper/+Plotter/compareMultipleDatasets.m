function compareMultipleDatasets(scanValsCell, meanValsCell, stderrValsCell, varargin)
%% compareMultipleDatasets
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Compares multiple datasets with error bars.
%
% Inputs:
%   scanValsCell   - cell array of x-values for each dataset
%   meanValsCell   - cell array of mean y-values for each dataset
%   stderrValsCell - cell array of std/error values for each dataset
%
% Name-Value Pair Arguments:
%   'FigNum', 'FontName', 'MarkerSize', 'LineWidth', 'CapSize',
%   'YLim', 'Labels', 'Title', 'XLabel', 'YLabel',
%   'SkipSaveFigures', 'SaveFileName', 'SaveDirectory'
%
% Notes:
%   Optional notes, references.

    % --- Parse inputs ---
    p = inputParser;
    addParameter(p, 'FigNum', 1, @isnumeric);
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'MarkerSize', 6, @isnumeric);
    addParameter(p, 'LineWidth', 1.5, @isnumeric);
    addParameter(p, 'CapSize', 5, @isnumeric);
    addParameter(p, 'YLim', [], @isnumeric);
    addParameter(p, 'Labels', {}, @iscell);
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'XLabel', '', @ischar);
    addParameter(p, 'YLabel', '', @ischar);
    addParameter(p, 'SkipSaveFigures', true, @islogical);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'SaveFileName', 'figure.fig', @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Default labels ---
    nDatasets = numel(scanValsCell);
    if isempty(opts.Labels)
        opts.Labels = arrayfun(@(i) sprintf('Dataset %d',i), 1:nDatasets, 'UniformOutput', false);
    end

    % --- Marker/line style cycle ---
    markerList = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+'};
    lineList   = {'-', '--', ':', '-.'};

    % --- Plot ---
    fig = figure(opts.FigNum); clf;
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);
    hold on;

    for i = 1:nDatasets
        marker = markerList{mod(i-1, numel(markerList)) + 1};
        lineStyle = lineList{mod(i-1, numel(lineList)) + 1};
        styleStr = [marker lineStyle];

        if isempty(stderrValsCell{i})
            plot(scanValsCell{i}, meanValsCell{i}, styleStr, ...
                'MarkerSize', opts.MarkerSize, 'LineWidth', opts.LineWidth, ...
                'DisplayName', opts.Labels{i});
        else
            errorbar(scanValsCell{i}, meanValsCell{i}, stderrValsCell{i}, styleStr, ...
                'MarkerSize', opts.MarkerSize, 'LineWidth', opts.LineWidth, 'CapSize', opts.CapSize, ...
                'DisplayName', opts.Labels{i});
        end
    end

    hold off;
    ax = gca;
    axisFontSize = 14;
    titleFontSize = 16;
    set(ax, 'FontName', opts.FontName, 'FontSize', axisFontSize);

    if ~isempty(opts.YLim)
        ylim(opts.YLim);
    end

    xlabel(opts.XLabel, 'Interpreter', 'latex', 'FontSize', axisFontSize);
    ylabel(opts.YLabel, 'Interpreter', 'latex', 'FontSize', axisFontSize);
    title(opts.Title, 'FontName', opts.FontName, 'FontSize', titleFontSize);
    legend('Location', 'best'); 
    grid on;

    % --- Save figure ---
    Plotter.saveFigure(fig, ...
        'SaveFileName', opts.SaveFileName, ...
        'SaveDirectory', opts.SaveDirectory, ...
        'SkipSaveFigures', opts.SkipSaveFigures);
end