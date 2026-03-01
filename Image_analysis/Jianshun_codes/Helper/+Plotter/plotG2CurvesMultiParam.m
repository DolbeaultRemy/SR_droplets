function plotG2CurvesMultiParam(results, scan_parameter_values, scan_reference_values, param1ValuesToPlot, param2ValuesToPlot, varargin)
%% plotG2CurvesMultiParam
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plot g²(θ) curves for selected param1 and param2 values in a phase diagram
%
% Notes:
%   Optional notes, references.

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', 'g² Curves', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Param1Name', 'Param1', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Param2Name', 'Param2', @(x) ischar(x) || isstring(x));
    addParameter(p, 'HighlightEvery', 10, @(x) isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'G2CurvesPhaseDiagram.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Convert reference/parameter values to numeric arrays ---
    scanRefArray   = cell2mat(scan_reference_values(:));   % nParams × 2
    scanParamArray = cell2mat(scan_parameter_values(:));   % nTotal × 2

    % --- Select subset based on param1ValuesToPlot and param2ValuesToPlot ---
    keepIdx = ismembertol(scanRefArray(:,1), param1ValuesToPlot, 1e-6) & ...
              ismembertol(scanRefArray(:,2), param2ValuesToPlot, 1e-6);
    scan_reference_values_plot = scan_reference_values(keepIdx);
    nParamsPlot = numel(scan_reference_values_plot);

    % --- Warn about missing values ---
    existingP1 = param1ValuesToPlot(ismembertol(param1ValuesToPlot, unique(scanRefArray(:,1)), 1e-6));
    missingP1  = setdiff(param1ValuesToPlot, existingP1);
    existingP2 = param2ValuesToPlot(ismembertol(param2ValuesToPlot, unique(scanRefArray(:,2)), 1e-6));
    missingP2  = setdiff(param2ValuesToPlot, existingP2);
    if ~isempty(missingP1)
        warning('The following %s values were not found: %s', opts.Param1Name, num2str(missingP1));
    end
    if ~isempty(missingP2)
        warning('The following %s values were not found: %s', opts.Param2Name, num2str(missingP2));
    end

    param1Plot = sort(existingP1, 'descend');
    param2Plot = sort(existingP2);

    % --- Extract theta values ---
    theta       = results.theta_values / pi;
    Nhighlight  = opts.HighlightEvery;

    % --- Create figure ---
    if isempty(opts.FigNum)
        fig = figure;
    else
        fig = figure(opts.FigNum);
    end
    clf(fig);
    set(fig,'Color','w','Position',[100 100 1200 800]);
    t = tiledlayout(numel(param1Plot), numel(param2Plot), 'TileSpacing','compact','Padding','compact');
    title(t, opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);

    % --- Loop over selected param1/param2 values ---
    for r = 1:numel(param1Plot)
        P1 = param1Plot(r);
        for c = 1:numel(param2Plot)
            P2 = param2Plot(c);

            ax = nexttile((r-1)*numel(param2Plot) + c); hold(ax,'on');

            % Find index of this (param1, param2)
            k = find(cellfun(@(v) all(abs(v - [P1 P2]) < 1e-6), scan_reference_values_plot), 1);
            if isempty(k), continue; end

            % Extract g² curves for this pair
            G = results.g2_curves{k}; % [N_reps × Nθ]

            % --- Plot all repetitions in light grey ---
            plot(ax, theta, G', 'Color', [0.7 0.7 0.7, 0.5]);

            % --- Highlight every Nth repetition ---
            idx = Nhighlight:Nhighlight:size(G,1);
            for j = idx
                plot(ax, theta, G(j,:), 'Color', [0.3 0.3 0.3, 1], 'LineWidth', 1.5);
            end

            % --- Mean + SEM shading ---
            mu = mean(G,1,'omitnan');
            se = std(G,0,1,'omitnan')/sqrt(size(G,1));
            fill(ax, [theta fliplr(theta)], [mu-se fliplr(mu+se)], [0.2 0.4 0.8], ...
                'FaceAlpha',0.2,'EdgeColor','none');

            % --- Mean curve ---
            plot(ax, theta, mu, 'b-', 'LineWidth', 2);

            % --- Vertical reference lines at pi/3 and 2pi/3 ---
            xlines = [1/3 2/3];
            for xl = xlines
                xline(ax, xl, 'k--', 'LineWidth', 1.5, 'Alpha', 0.5);
            end

            % --- Axes formatting ---
            grid(ax,'on');
            xlabel(ax, '\delta\theta / \pi', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
            ylabel(ax, 'g^2(\theta)', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
            title(ax, sprintf('%s=%.3g, %s=%.1f', opts.Param1Name, P1, opts.Param2Name, P2), ...
                  'FontName', opts.FontName, 'FontSize', opts.FontSize);
            set(ax, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        end
    end

    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end
