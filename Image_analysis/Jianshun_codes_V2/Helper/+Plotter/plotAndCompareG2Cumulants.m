function plotAndCompareG2Cumulants(results1, results2, varargin)
%% plotAndCompareG2Cumulants
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Compare first four cumulants of g²(θ) for two datasets.
%
% Inputs:
%   results1, results2 : g2_analysis_results objects
%
% Notes:
%   Dataset 1 = Solid line
%   Dataset 2 = Dashed line
%   Marker shape encodes θ value (consistent across datasets)
%   Dataset colors are distinct

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', 'g² Cumulants Comparison', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', 'Scan Parameter', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'G2CumulantsComparison.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'DesiredTheta', [pi/5 pi/3 pi/2], @(x) isnumeric(x) && isvector(x)); % <--- new line
    parse(p, varargin{:});
    opts = p.Results;


    % --- Helper to extract cumulant array ---
    function [xvals, kappa, thIdx] = extractData(results)
        N_params = numel(results.g2_cumulants);
        N_theta  = size(results.g2_cumulants{1},1);
        xvals    = results.summaries.scanParamVal;
        kappa    = zeros(N_theta, N_params, 4);
        for i = 1:N_params
            kappa(:,i,:) = results.g2_cumulants{i};
        end
        thetaVals   = results.theta_values;
        desiredTheta = opts.DesiredTheta;   % use user-specified or default
        [~, thIdx]  = arrayfun(@(t) min(abs(thetaVals - t)), desiredTheta);
    end

    % --- Extract data for both results ---
    [x1, kappa1, thIdx] = extractData(results1);
    [x2, kappa2, ~]     = extractData(results2);

    % --- Legend labels for θ ---
    thetaLabels = arrayfun(@(t) sprintf('\\pi/%g', round(pi/t)), opts.DesiredTheta, 'UniformOutput', false);

    % --- Line styles for datasets ---
    lineStyles = {'-', '--'};  % solid for dataset1, dashed for dataset2

    % --- Colors per θ ---
    %{
    cmapFull = Colormaps.coolwarm(256);
    Ntheta = numel(thIdx);
    startColor = cmapFull(1,:); 
    endColor   = cmapFull(end,:);
    thetaColors = [linspace(startColor(1), endColor(1), Ntheta)', ...
                   linspace(startColor(2), endColor(2), Ntheta)', ...
                   linspace(startColor(3), endColor(3), Ntheta)'];
    %}
    cmapFull = Colormaps.coolwarm(256);
    Ntheta   = numel(thIdx);
    % Split into two halves: blue side (1:128), red side (129:256)
    blueSide = cmapFull(1:100, :);    % cut before gray
    redSide  = cmapFull(157:end, :);  % cut after gray
    % Concatenate to avoid center
    cmapTrimmed = [blueSide; redSide];
    % Now sample evenly from this trimmed colormap
    indices = round(linspace(1, size(cmapTrimmed,1), Ntheta));
    thetaColors = cmapTrimmed(indices, :);


    % --- Create figure ---
    if isempty(opts.FigNum), fig = figure; else, fig = figure(opts.FigNum); end
    clf(fig);
    set(fig,'Color','w','Position',[100 100 950 750]);
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    sgtitle(opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize+4);

    cumulLabels = {'\kappa_1','\kappa_2','\kappa_3','\kappa_4'};
    cumulTitles = {'Mean','Variance','Skewness','Binder Cumulant'};
    
    for k = 1:4
        ax = nexttile; hold(ax,'on');
    
        for idx = 1:Ntheta
            th = thIdx(idx);
            c = thetaColors(idx,:);
            
            % Dataset 1 -> solid line + marker
            plot(ax, x1, squeeze(kappa1(th,:,k)), ...
                'LineStyle', lineStyles{1}, 'Color', c, ...
                'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 7, ... % markerList = {'o','s','^','d','v','>','<','p','h'}; % more markers if needed
                'MarkerFaceColor', c);

            % Dataset 2 -> dashed line + marker
            plot(ax, x2, squeeze(kappa2(th,:,k)), ...
                'LineStyle', lineStyles{2}, 'Color', c, ...
                'LineWidth', 2, 'Marker', 's', 'MarkerSize', 7, ...
                'MarkerFaceColor', c);
        end
    
        ylabel(ax, cumulLabels{k}, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        xlabel(ax, opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        title(ax, cumulTitles{k}, 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);
        grid(ax,'on'); set(ax,'FontName',opts.FontName,'FontSize',opts.FontSize);

        % --- Line-only legend for θ ---
        thetaHandles = gobjects(Ntheta,1);
        for m = 1:Ntheta
            thetaHandles(m) = plot(ax, NaN, NaN, ...
                'LineStyle','-', 'Color',thetaColors(m,:), 'LineWidth',6, 'Marker','none');
        end
        legend(ax, thetaHandles, thetaLabels, 'Location','northeast', 'FontSize', opts.FontSize-2);
    end

    % --- Dataset legend (solid/dashed + marker) ---
    axLegend = axes(fig,'Position',[0 0.91 1 0.05],'Visible','off'); 
    hold(axLegend,'on');
    
    % Dummy plots with line style + marker in black
    d1 = plot(axLegend, NaN, NaN, 'LineStyle', lineStyles{1}, ...
        'Color','k', 'Marker','o', 'MarkerFaceColor','k', 'MarkerSize',7, 'LineWidth',2);
    d2 = plot(axLegend, NaN, NaN, 'LineStyle', lineStyles{2}, ...
        'Color','k', 'Marker','s', 'MarkerFaceColor','k', 'MarkerSize',7, 'LineWidth',2);
    
    legend(axLegend,[d1 d2],{'D -> S','S -> D'}, ...
        'FontName', opts.FontName, 'FontSize', opts.FontSize-2, ...
        'Orientation','horizontal', 'Box','off', 'Location','north');

    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end
