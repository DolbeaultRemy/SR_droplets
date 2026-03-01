function plotSinglePCAResults(pcaResults, scan_parameter_values, scan_reference_values,  varargin)
%% plotSinglePCAResults
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots PCA results.
%
% Inputs:
%   pcaResults - struct returned by computePCAfromImages
%   scan_parameter_values, scan_reference_values
%   varargin   - name-value pairs (same as plotG2 plus 'FigNumRange','MaxPCToPlot')
%
% Notes:
%   Optional notes, references.

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'XLabel', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'Colormap', @Colormaps.coolwarm);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'FigNumRange', [], @(x) isnumeric(x) && all(x>0));
    parse(p, varargin{:});
    opts = p.Results;

    Nx         = pcaResults.Nx;
    Ny         = pcaResults.Ny;
    coeff      = pcaResults.coeff;
    score      = pcaResults.score;
    explained  = pcaResults.explained;

    raw_scan_param_vals    = scan_parameter_values;
    unique_scan_param_vals = scan_reference_values;
    numGroups              = numel(unique_scan_param_vals);
    colors                 = lines(numGroups);

    % --- Figure numbering setup ---
    if isempty(opts.FigNumRange)
        figCount = 1;
        figNums = [];
    else
        figNums = opts.FigNumRange;
        figCount = 1;
    end

    figPos = [100 100 950 750];

    %% --- Figure 1: PC1 Image ---
    pc1_image = reshape(coeff(:,1), Nx, Ny);
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position', figPos);
    imagesc(pc1_image); axis image off; colormap(opts.Colormap()); colorbar;
    title(sprintf('First Principal Component (PC1) Image - Explains %.2f%% Variance', explained(1)), ...
        'FontName', opts.FontName, 'FontSize', opts.FontSize + 2);
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', 'PC1_Image.fig', 'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;

    %% --- Figure 2: PC1 Scores Distribution Scatterplot ---
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position', figPos); hold on;
    for g = 1:numGroups
        idx = raw_scan_param_vals == unique_scan_param_vals(g);
        scatter(repmat(unique_scan_param_vals(g), sum(idx),1), score(idx,1), 36, colors(g,:), 'filled');
    end
    xlabel(opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel('PC1 Score', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title('Evolution of PC1 Scores', 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2);
    grid on;
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', 'PC1_Scatter.fig', 'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;

    %% --- Figure 3: PC1 Scores Distribution Histograms ---
    numTiles = min(6, numGroups);  % show up to 6 groups
    tileIndices = round(linspace(1, numGroups, numTiles));

    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position', figPos);
    tLayout = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    for t = 1:numTiles
        g = tileIndices(t);
        idx = raw_scan_param_vals == unique_scan_param_vals(g);
        data = score(idx,1);

        nexttile;
        histogram(data, 'Normalization', 'pdf', 'FaceColor', colors(g,:), 'FaceAlpha', 0.3);
        hold on;
        [f, xi] = ksdensity(data);
        plot(xi, f, 'Color', colors(g,:), 'LineWidth', 2);
        yl = ylim;
        plot([median(data) median(data)], yl, 'k--', 'LineWidth', 1);
        xlabel('PC1 Score', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        ylabel('Probability', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        title(sprintf('Control = %g', unique_scan_param_vals(g)), 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        grid on;
    end
    sgtitle('PC1 Score Distributions', 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2);
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', 'PC1_Distributions.fig', 'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;

    %% --- Figure 4: PC1 Scores Distribution Boxplot ---
    % Construct group labels explicitly
    groupLabels = arrayfun(@num2str, raw_scan_param_vals, 'UniformOutput', false);
    
    % Create categorical variable with specified order
    groupCats = categorical(groupLabels, ...
                            arrayfun(@num2str, unique_scan_param_vals, 'UniformOutput', false), ...
                            'Ordinal', true);
    
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position', figPos);
    
    % Plot boxplot with categorical groups
    boxplot(score(:,1), groupCats);
    
    xlabel(opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel('PC1 Score', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title('PC1 Score Boxplots by Group', 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2);
    grid on;
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', 'PC1_Boxplot.fig', 'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;


    %% --- Figure 5: PC1 Scores Distribution Mean ± SEM ---
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position', figPos);
    meanScores = arrayfun(@(g) mean(score(raw_scan_param_vals == g,1)), unique_scan_param_vals);
    semScores  = arrayfun(@(g) std(score(raw_scan_param_vals == g,1))/sqrt(sum(raw_scan_param_vals == g)), unique_scan_param_vals);
    errorbar(unique_scan_param_vals, meanScores, semScores, '--o', 'LineWidth', 2);
    xlabel(opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel('Mean PC1 Score ± SEM', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title('Mean ± SEM of PC1 Scores', 'FontName', opts.FontName, 'FontSize', opts.FontSize + 2);
    grid on;
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', 'PC1_MeanSEM.fig', 'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;

    %% --- Figure 6: PC1 Scores Distribution Cumulants ---
    kappas = cell2mat(arrayfun(@(g) {Calculator.computeCumulants(score(raw_scan_param_vals == g,1), 4)}, unique_scan_param_vals));
    
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig,'Color','w','Position',[100 100 950 750]);
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    title(t, 'Cumulants of PC1 Scores', ...
        'FontName', opts.FontName, 'FontSize', opts.FontSize+4);

    cumulLabels = {'\kappa_1','\kappa_2','\kappa_3','\kappa_4'};
    cumulTitles = {'Mean','Variance','Skewness','Binder Cumulant'};

    for k = 1:4
        ax = nexttile; hold(ax,'on');
        plot(ax, unique_scan_param_vals, kappas(:, k), '-o', ...
            'Color', [0.2 0.4 0.7], 'LineWidth', 2, 'MarkerSize', 8, ...
            'MarkerFaceColor', [0.2 0.4 0.7]);
        ylabel(ax, cumulLabels{k}, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        xlabel(ax, opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        title(ax, cumulTitles{k}, 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);
        grid(ax,'on'); set(ax,'FontName',opts.FontName,'FontSize',opts.FontSize);
    end
    
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', 'PC1_BinderCumulant.fig', 'SaveDirectory', opts.SaveDirectory);
    end
end
