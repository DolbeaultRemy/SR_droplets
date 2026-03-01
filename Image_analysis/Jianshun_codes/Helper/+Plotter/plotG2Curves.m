function plotG2Curves(results, varargin)
%% plotG2Curves
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plot raw g²(θ) curves with mean, SEM, and highlights
%
% Inputs:
%   results - struct with fields:
%       g2_curves            - cell array of [N_reps × Nθ] matrices
%       theta_values         - vector of theta values (radians)
%       g2_mean              - [N_params × Nθ] mean
%       g2_error             - [N_params × Nθ] SEM
%       scan_parameter_values- vector of scan parameters
%
% Notes:
%   Optional notes, references. 

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', 'g² Curves', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', '\theta / \pi', @(x) ischar(x) || isstring(x));
    addParameter(p, 'YLabel', 'g^2(\theta)', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'G2Curves.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'TileTitlePrefix', 'Control Parameter', @(x) ischar(x) || isstring(x));
    addParameter(p, 'TileTitleSuffix', '', @(x) ischar(x) || isstring(x)); % NEW
    addParameter(p, 'HighlightEvery', 10, @(x) isnumeric(x) && isscalar(x) && x>=1);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Extract data ---
    N_params   = numel(results.g2_curves);
    theta      = results.theta_values / pi;
    g2_mean    = results.g2_mean;
    g2_err     = results.g2_error;
    Nhighlight = opts.HighlightEvery;

    % --- Create figure ---
    if isempty(opts.FigNum), fig = figure; else, fig = figure(opts.FigNum); end
    clf(fig);
    set(fig,'Color','w','Position',[100 100 950 750]);
    t = tiledlayout('TileSpacing','compact','Padding','compact');
    title(t, opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);

    % --- Loop over scan parameters ---
    for i = 1:N_params
        ax = nexttile; hold(ax,'on');

        G = results.g2_curves{i}; % [N_reps × Nθ]

        % --- Plot all repetitions in light grey ---
        plot(ax, theta, G', 'Color', [0.7 0.7 0.7, 0.5]);

        % --- Highlight every Nth curve ---
        idx = Nhighlight:Nhighlight:size(G,1);
        for j = idx
            plot(ax, theta, G(j,:), 'Color', [0.3 0.3 0.3, 1], 'LineWidth', 1.5);
        end

        % --- Mean + SEM shading ---
        mu = g2_mean{i};
        se = g2_err{i};
        fill(ax, [theta fliplr(theta)], [mu-se fliplr(mu+se)], ...
            [0.2 0.4 0.8], 'FaceAlpha',0.2, 'EdgeColor','none');

        % --- Mean curve ---
        plot(ax, theta, mu, 'b-', 'LineWidth', 2);

        % --- Vertical reference lines at pi/3 and 2pi/3 ---
        xlines = [1/3 2/3];
        for xl = xlines
            xline(ax, xl, 'k--', 'LineWidth', 1.5, 'Alpha', 0.5);
        end

        % --- Axes formatting ---
        grid(ax,'on');
        xlabel(ax, opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        ylabel(ax, opts.YLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        title(ax, sprintf('%s=%.3g%s', opts.TileTitlePrefix, results.scan_parameter_values(i), opts.TileTitleSuffix), ...
                'FontName', opts.FontName, 'FontSize', opts.FontSize);
        set(ax, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    end

    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end