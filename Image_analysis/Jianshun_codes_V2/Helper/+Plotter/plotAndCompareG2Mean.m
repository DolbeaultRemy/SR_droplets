function plotAndCompareG2Mean(results1, results2, varargin)
%% plotAndCompareG2Mean
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Compare first cumulant (mean) of g²(θ) for two datasets.
%   Plots mean ± standard deviation of g² values at the desired θ
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
    addParameter(p, 'Title', 'Mean Comparison', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', 'Scan Parameter', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLabel', 'Mean (\kappa_1)', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLimits', [], @(x) isnumeric(x) && numel(x)==2);
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'G2MeanComparison.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'DesiredTheta', [pi/12 pi/6 pi/3 pi/2 2*pi/3 5*pi/6], @(x) isnumeric(x) && isvector(x));
    parse(p, varargin{:});
    opts = p.Results;

    % --- Helper to extract mean, STD, and SEM at desired θ ---
    function [xvals, yMean, yStd, ySEM, thIdx] = extractData(results)
        xvals        = results.scan_parameter_values(:); % [N_params × 1]
        N_params     = numel(xvals);
        thetaVals    = results.theta_values;             % all θ
        desiredTheta = opts.DesiredTheta;
        [~, thIdx]   = arrayfun(@(t) min(abs(thetaVals - t)), desiredTheta);
        Ndesired     = numel(desiredTheta);
    
        yMean        = zeros(N_params, Ndesired); % [scan × θ]
        yStd         = zeros(N_params, Ndesired);
        ySEM         = zeros(N_params, Ndesired);
    
        for i = 1:N_params
            G = results.g2_curves{i};                         % [N_reps × N_theta] for this scan parameter
            for j = 1:Ndesired
                th = thIdx(j);
                data = G(:, th);                              % all reps for this θ at this scan param
                yMean(i,j) = mean(data);                      % mean across reps
                yStd(i,j)  = std(data,0,1);                   % std across reps
                ySEM(i,j)  = std(data,0,1)/sqrt(numel(data)); % SEM
            end
        end
    end


    % --- Extract data for both datasets ---
    [x1, yMean1, yStd1, ySEM1, thIdx] = extractData(results1);
    [x2, yMean2, yStd2, ySEM2, ~]     = extractData(results2);

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
    
    ax = axes(fig); hold(ax,'on'); grid(ax,'on');
    % Set axes font for labels and numbers
    set(ax, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    if ~isempty(opts.YLimits)
        ylim(ax, opts.YLimits);
    end
    xlabel(ax, opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel(ax, opts.YLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    title(ax, opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);

    % --- Plot each θ with SEM ---
    for idx = 1:Ntheta
        c = thetaColors(idx,:);

        % Dataset 1 -> solid line + marker + SEM
        errorbar(ax, x1, yMean1(:,idx), ySEM1(:,idx), 'LineStyle', lineStyles{1}, ...
            'Color', c, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', c);

        % Dataset 2 -> dashed line + marker + SEM
        errorbar(ax, x2, yMean2(:,idx), ySEM2(:,idx), 'LineStyle', lineStyles{2}, ...
            'Color', c, 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 7, 'MarkerFaceColor', c);
    end
    
    % --- Line-only legend for θ ---
    thetaHandles = gobjects(Ntheta,1);
    for m = 1:Ntheta
        thetaHandles(m) = plot(ax, NaN, NaN, 'LineStyle','-', 'Color',thetaColors(m,:), ...
            'LineWidth',6, 'Marker','none');
    end
    legend(ax, thetaHandles, thetaLabels, 'Location','northeast', 'FontSize', opts.FontSize-2);

    % --- Dataset legend (solid/dashed + marker) ABOVE AXES ---
    axLegend = axes(fig,'Position',[0.65 0.65 0.25 0.15],'Visible','off'); 
    hold(axLegend,'on');
    
    % Create dummy plots for legend with line style + marker in black
    d1 = plot(axLegend, NaN, NaN, 'LineStyle', lineStyles{1}, ...
        'Color','k', 'Marker','o', 'MarkerFaceColor','k', 'MarkerSize',7, 'LineWidth',2);
    d2 = plot(axLegend, NaN, NaN, 'LineStyle', lineStyles{2}, ...
        'Color','k', 'Marker','s', 'MarkerFaceColor','k', 'MarkerSize',7, 'LineWidth',2);
    
    legend(axLegend,[d1 d2],{'D -> S','S -> D'}, ...
        'FontName', opts.FontName, 'FontSize', opts.FontSize+2, ...
        'Orientation','vertical', 'Box','off');



    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end
