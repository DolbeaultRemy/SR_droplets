function plotG2Cumulants(results, varargin)
%% plotG2Cumulants
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plot first four cumulants of g²(θ) vs scan parameter
%
% Inputs:
%   img        - original 2D image
%   patchProps - struct array from detectPatches
%   results    - struct array from extractAndClassifyShapes
%   params     - parameter struct
%   xStart, yStart - offsets for cropped regions
%   figTitle   - title string (e.g., 'BayesOpt Candidate')
%
% Notes:
%   Optional notes, references. 

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', 'g² Cumulants', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', 'Scan Parameter', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'G2Cumulants.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Extract cumulant data from full distribution ---
    N_params = numel(results.g2_cumulants);
    N_theta  = size(results.g2_cumulants{1},1);
    xvals    = results.summaries.scanParamVal;

    kappa = zeros(N_theta, N_params, 4); % θ × scan × cumulant
    for i = 1:N_params
        cumul_mat = results.g2_cumulants{i}; % [N_theta × 4]
        kappa(:,i,:) = cumul_mat;
    end

    % --- Select specific theta indices: 0, pi/6, pi/3, pi/2, 2*pi/3, pi ---
    thetaVals = results.theta_values; % assume radians
    desiredTheta = [pi/12 pi/6 pi/3 pi/2 2*pi/3 5*pi/6]; 
    [~, thIdx] = arrayfun(@(t) min(abs(thetaVals - t)), desiredTheta);

    % --- Legend labels ---
    thetaLabels = {'\pi/12','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6'};

    % --- Colormap (sky) ---
    cmap = Colormaps.coolwarm(numel(thIdx));

    % --- Create figure ---
    if isempty(opts.FigNum), fig = figure; else, fig = figure(opts.FigNum); end
    clf(fig);
    set(fig,'Color','w','Position',[100 100 950 750]);
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    title(t, opts.Title, 'FontName', opts.FontName, 'FontSize', opts.FontSize+4);

    cumulLabels = {'\kappa_1','\kappa_2','\kappa_3','\kappa_4'};
    cumulTitles = {'Mean','Variance','Skewness','Binder Cumulant'};

    for k = 1:4
        ax = nexttile; hold(ax,'on');
        for idx = 1:numel(thIdx)
            th = thIdx(idx);
            plot(ax, xvals, squeeze(kappa(th,:,k)), '-o', ...
                'Color', cmap(idx,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', cmap(idx,:));
        end
        ylabel(ax, cumulLabels{k}, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        xlabel(ax, opts.XLabel, 'FontName', opts.FontName, 'FontSize', opts.FontSize);
        title(ax, cumulTitles{k}, 'FontName', opts.FontName, 'FontSize', opts.FontSize+2);
        grid(ax,'on'); set(ax,'FontName',opts.FontName,'FontSize',opts.FontSize);

        % --- Add legend ---
        legend(ax, thetaLabels, 'Location', 'best', 'FontSize', opts.FontSize-2);

    end

    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end
