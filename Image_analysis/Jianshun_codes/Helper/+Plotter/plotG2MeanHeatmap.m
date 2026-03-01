function plotG2MeanHeatmap(results, theta_query, scan_reference_values, varargin)
%% plotG2MeanHeatmap
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plots a heatmap of mean g² values at a specified theta across the phase
%   diagram.
%
% Notes:
%   Optional notes, references.

    % --- Parse optional inputs ---
    p = inputParser;
    addParameter(p, 'FigNum', []);
    addParameter(p, 'Colormap', @jet);
    addParameter(p, 'CLim', []);
    addParameter(p, 'ColorScale', 'linear', @(x) any(validatestring(x,{'linear','log'})));
    addParameter(p, 'XLabel', 'Param1');
    addParameter(p, 'YLabel', 'Param2');
    addParameter(p, 'Title', 'Mean g²');
    addParameter(p, 'FontName', 'Arial');
    addParameter(p, 'FontSize', 14);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'g2MeanHeatmap.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Extract data at requested theta ---
    theta_vals      = results.theta_values;
    [~, idx_theta]  = min(abs(theta_vals - theta_query));  % closest index
    N_total         = numel(results.g2_curves);
    dataCell        = cell(N_total,1);
    for i = 1:N_total
        G           = results.g2_curves{i};       % [N_reps × Nθ]
        dataCell{i} = G(:, idx_theta);   % all repetitions at chosen theta
    end

    % --- Convert reference/parameter values to numeric arrays ---
    scanRefArray    = cell2mat(scan_reference_values(:));   % nParams × 2

    % --- Determine unique X and Y values ---
    XValues         = sort(unique(scanRefArray(:,2)), 'ascend'); % Param1 → X
    YValues         = sort(unique(scanRefArray(:,1)), 'ascend'); % Param2 → Y
    nX              = numel(XValues);
    nY              = numel(YValues);

    % --- Preallocate data matrix ---
    data_matrix     = NaN(nY, nX);

    % --- Fill matrix by looping over parameter pairs ---
    for k = 1:N_total
        X                    = scanRefArray(k,2);
        Y                    = scanRefArray(k,1);
    
        row                  = find(ismembertol(YValues, Y, 1e-6));
        col                  = find(ismembertol(XValues, X, 1e-6));
        
        if isempty(row) || isempty(col), continue; end
    
        % Get g² values for this (Y,X)
        G                    = results.g2_curves{k};      % [N_reps × Nθ]
        vals                 = G(:, idx_theta);           % repetitions at chosen θ
    
        % Store mean
        data_matrix(row,col) = mean(vals, 'omitnan');
    end

    % --- Create figure ---
    if isempty(opts.FigNum)
        fig = figure;
    else
        fig = figure(opts.FigNum);
    end
    clf(fig);
    set(fig, 'Color', 'w', 'Position', [100 100 950 750]);

    % --- Plot heatmap ---
    imagesc(XValues, YValues, data_matrix);
    set(gca, ...
        'FontName', opts.FontName, ...
        'FontSize', opts.FontSize, ...
        'YDir', 'normal', ...
        'ColorScale', opts.ColorScale);

    % --- Labels and title ---
    xlabel(opts.XLabel, 'Interpreter', 'tex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);
    ylabel(opts.YLabel, 'Interpreter', 'tex', 'FontName', opts.FontName, 'FontSize', opts.FontSize);

    % --- Compute fraction of pi for theta label ---
    theta_factor = theta_query / pi;  
    if abs(round(1/theta_factor) - 1/theta_factor) < 1e-6  
        theta_str = sprintf('\\pi/%d', round(1/theta_factor));
    elseif abs(theta_factor - round(theta_factor)) < 1e-6   
        theta_str = sprintf('%d\\pi', round(theta_factor));
    else
        theta_str = sprintf('%.3g \\pi', theta_factor);  
    end

    % --- Set title ---
    title(sprintf('%s | Mean g^{(2)}(\\delta\\theta = %s)', ...
          opts.Title, theta_str), ...
          'FontName', opts.FontName, 'FontSize', opts.FontSize + 2, 'FontWeight', 'bold');
    
    % --- Colorbar ---
    colormap(feval(opts.Colormap));
    if ~isempty(opts.CLim)
        caxis(opts.CLim);
    end
    colorbar;

    % --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end