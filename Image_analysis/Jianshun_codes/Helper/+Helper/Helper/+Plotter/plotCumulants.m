function plotCumulants(scan_vals, cumulant_data, varargin)
%% plotCumulants
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots the first four cumulants vs. a scan parameter
%
% Inputs:
%   scan_vals      - array of scan parameter values
%   cumulant_data  - cell array of cumulants:
%
% Notes:
%   Optional notes, references.

    % --- Parse optional name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'XLabel', 'Scan Parameter', @ischar);
    addParameter(p, 'FigNum', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'MarkerSize', 8, @isnumeric);
    addParameter(p, 'LineWidth', 2, @isnumeric);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'cumulants.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Extract cumulant data ---
    mean_vals         = cumulant_data{1};
    var_vals          = cumulant_data{2};
    skew_vals         = cumulant_data{3};
    fourth_order_vals = cumulant_data{4};

    % --- Figure setup ---
    fig = figure(opts.FigNum); clf;
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);

    axisFontSize  = 14;
    labelFontSize = 16;
    titleFontSize = 16;

    tLayout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % Define style
    plotColor = [0.2 0.4 0.7];

    % --- Mean ---
    nexttile; hold on;
    plot(scan_vals, mean_vals, '-o', ...
        'Color', plotColor, 'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, ...
        'MarkerFaceColor', plotColor);
    % add error bars on top of styled plot
    plot(scan_vals, mean_vals, '-o', ...
            'Color', plotColor, 'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, ...
            'MarkerFaceColor', plotColor);
    title('Mean', 'FontSize', titleFontSize, 'FontWeight', 'bold');
    xlabel(opts.XLabel, 'FontSize', labelFontSize); 
    ylabel('\kappa_1', 'FontSize', labelFontSize);
    set(gca, 'FontSize', axisFontSize, 'FontName', opts.FontName);
    grid on;

    % --- Variance ---
    nexttile; hold on;
    plot(scan_vals, var_vals, '-o', ...
        'Color', plotColor, 'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, ...
        'MarkerFaceColor', plotColor);
    title('Variance', 'FontSize', titleFontSize, 'FontWeight', 'bold');
    xlabel(opts.XLabel, 'FontSize', labelFontSize); 
    ylabel('\kappa_2', 'FontSize', labelFontSize);
    set(gca, 'FontSize', axisFontSize, 'FontName', opts.FontName);
    grid on;

    % --- Skewness ---
    nexttile; hold on;
    plot(scan_vals, skew_vals, '-o', ...
        'Color', plotColor, 'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, ...
        'MarkerFaceColor', plotColor);
    title('Skewness', 'FontSize', titleFontSize, 'FontWeight', 'bold');
    xlabel(opts.XLabel, 'FontSize', labelFontSize); 
    ylabel('\kappa_3', 'FontSize', labelFontSize);
    set(gca, 'FontSize', axisFontSize, 'FontName', opts.FontName);
    grid on;

    % --- Binder Cumulant ---
    nexttile; hold on;
    plot(scan_vals, fourth_order_vals, '-o', ...
        'Color', plotColor, 'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, ...
        'MarkerFaceColor', plotColor);
    title('Binder Cumulant', 'FontSize', titleFontSize, 'FontWeight', 'bold');
    xlabel(opts.XLabel, 'FontSize', labelFontSize); 
    ylabel('\kappa_4', 'FontSize', labelFontSize);
    set(gca, 'FontSize', axisFontSize, 'FontName', opts.FontName);
    grid on;

    % --- Super title ---
    if ~isempty(opts.Title)
        sgtitle(opts.Title, 'FontName', opts.FontName, 'FontWeight', 'bold', 'FontSize', titleFontSize);
    end

    % --- Save figure ---
    Plotter.saveFigure(fig, ...
        'SaveFileName', opts.SaveFileName, ...
        'SaveDirectory', opts.SaveDirectory, ...
        'SkipSaveFigures', opts.SkipSaveFigures);

end
