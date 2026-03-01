function plotMultiplePCAResults(pcaResults, scan_parameter_values, scan_reference_values, varargin)
%% plotMultiplePCAResults
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots PCA results for multiple PCs.
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
    addParameter(p, 'MaxPCToPlot', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
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

    % --- Figure numbering setup ---
    if isempty(opts.FigNumRange)
        figCount = 1;
        figNums = [];
    else
        figNums = opts.FigNumRange;
        figCount = 1;
    end

    figPos = [100 100 950 750];

    %% --- Precompute score norms and group by scan parameter ---
    numPCsToPlot    = min(opts.MaxPCToPlot, size(coeff,2));
    scoreNorms      = sqrt(sum(score(:,1:numPCsToPlot).^2, 2));

    groupedNorms    = arrayfun(@(g) scoreNorms(raw_scan_param_vals == g), unique_scan_param_vals, 'UniformOutput', false);

    % Mean and SEM per group
    meanNorms       = cellfun(@mean, groupedNorms);
    semNorms        = cellfun(@(x) std(x)/sqrt(numel(x)), groupedNorms);

    % Cumulants (up to 4th order) per group
    kappas          = cell2mat(cellfun(@(x) Calculator.computeCumulants(x(:),4), groupedNorms, 'UniformOutput', false));    


    %% --- Figure 1: PC images ---
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position',  figPos);
    nRows = ceil(sqrt(numPCsToPlot));
    nCols = ceil(numPCsToPlot/nRows);
    tLayout = tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');

    for pc = 1:numPCsToPlot
        nexttile;
        pc_image = reshape(coeff(:,pc), Nx, Ny);
        imagesc(pc_image); axis image off;
        colormap(opts.Colormap());
        title(sprintf('PC%d (%.2f%%)', pc, explained(pc)), ...
            'FontName', opts.FontName, 'FontSize', opts.FontSize);
    end
    sgtitle(sprintf('Principal Component Images (1-%d)', numPCsToPlot), ...
        'FontName', opts.FontName, 'FontSize', opts.FontSize+2);
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);

    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', sprintf('PC1to%d_Images.fig', numPCsToPlot), ...
                           'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;

    %% --- Figure 2: Mean ± SEM of score norms ---
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig, 'Color', 'w', 'Position', figPos);

    errorbar(unique_scan_param_vals, meanNorms, semNorms, '--o', 'LineWidth', 2);
    
    xlabel(opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel(sprintf('‖Scores(1:%d)‖ ± SEM', numPCsToPlot), 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title('Mean ± SEM of Score Norms', 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);
    set(gca, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    grid on;

    if ~opts.SkipSaveFigures
        Plotter.saveFigure(fig, 'SaveFileName', sprintf('PC1to%d_Norm_MeanSEM.fig', numPCsToPlot), ...
                           'SaveDirectory', opts.SaveDirectory);
    end
    figCount = figCount + 1;

    %% --- Figure 3: Cumulants of score norms ---
    if ~isempty(figNums)
        fig = figure(figNums(figCount)); clf;
    else
        fig = figure; clf;
    end
    set(fig,'Color','w','Position',[100 100 950 750]);
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    title(t, sprintf('Cumulants of Score Norms (1-%d PCs)', numPCsToPlot), ...
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
        Plotter.saveFigure(fig, 'SaveFileName', sprintf('PC1to%d_Norm_Cumulants.fig', numPCsToPlot), ...
                           'SaveDirectory', opts.SaveDirectory);
    end

end
