function plotAverageSpectra(scan_parameter_values, spectral_analysis_results, varargin)
%% plotAverageSpectra
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plot averaged power, radial, and angular spectra for a scan
%
% Inputs:
%   scan_parameter_values      - array of scan parameter values
%   spectral_analysis_results  - struct with fields:
%       kx, ky, PS_all, k_rho_vals, S_k_all, theta_vals, S_theta_all
%
% Name-Value Pair Arguments:
%   'ScanParameterName', 'FigNum', 'ColormapPS', 'Font', 
%   'SaveFileName', 'SaveDirectory', 'SkipSaveFigures'
%
% Notes:
%   Dataset 1 = Solid line
%   Dataset 2 = Dashed line
%   Marker shape encodes θ value (consistent across datasets)
%   Dataset colors are distinct

    % --- Extract data from struct ---
    kx           = spectral_analysis_results.kx;
    ky           = spectral_analysis_results.ky;
    ps_list      = spectral_analysis_results.PS_all;
    k_rho_vals   = spectral_analysis_results.k_rho_vals;
    s_k_list     = spectral_analysis_results.S_k_all;
    theta_vals   = spectral_analysis_results.theta_vals;
    s_theta_list = spectral_analysis_results.S_theta_all;

    % --- Parse optional parameters ---
    p = inputParser;
    addParameter(p, 'ScanParameterName', 'ScanParameter', @ischar);
    addParameter(p, 'FigNum', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ColormapPS', Colormaps.coolwarm(), @(x) isnumeric(x) || ismatrix(x));
    addParameter(p, 'Font', 'Arial', @ischar);
    addParameter(p, 'SaveFileName', 'avgspectra.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    % --- Unique scan parameters ---
    [uniqueParams, ~, idx] = unique(scan_parameter_values);
    nParams = numel(uniqueParams);

    % --- Loop over each unique parameter ---
    for pIdx = 1:nParams
        currentParam = uniqueParams(pIdx);
        shotIndices  = find(idx == pIdx);
        nShots       = numel(shotIndices);

        % --- Initialize accumulators ---
        avgPS      = 0;
        avgS_k     = 0;
        avgS_theta = 0;

        % --- Sum over shots ---
        for j = 1:nShots
            avgPS      = avgPS      + ps_list{shotIndices(j)};
            avgS_k     = avgS_k     + s_k_list{shotIndices(j)};
            avgS_theta = avgS_theta + s_theta_list{shotIndices(j)};
        end

        % --- Average ---
        avgPS      = avgPS / nShots;
        avgS_k     = avgS_k / nShots;
        avgS_theta = avgS_theta / nShots;

        % ==== Plot ====
        fig = figure(opts.FigNum); clf;
        set(fig, 'Color', 'w', 'Position', [400 200 1200 400]);
        tLayout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

        axisFontSize = 14;
        titleFontSize = 16;

        % --- 1. Power Spectrum ---
        nexttile;
        imagesc(kx, ky, log(1 + avgPS));
        axis image;
        set(gca, 'FontSize', axisFontSize, 'YDir', 'normal');
        xlabel('k_x [\mum^{-1}]','Interpreter','tex','FontSize',axisFontSize,'FontName',opts.Font);
        ylabel('k_y [\mum^{-1}]','Interpreter','tex','FontSize',axisFontSize,'FontName',opts.Font);
        title('Average Power Spectrum','FontSize',titleFontSize,'FontWeight','bold');
        colormap(opts.ColormapPS);
        colorbar;

        % --- Annotate scan parameter ---
        if strcmp(opts.ScanParameterName,'ps_rot_mag_fin_pol_angle')
            txt = sprintf('%.1f^\\circ', currentParam);
        else
            txt = sprintf('%.2f G', currentParam);
        end
        text(0.975,0.975,txt,'Color','white','FontWeight','bold','FontSize',axisFontSize, ...
            'Interpreter','tex','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');

        % --- 2. Radial Spectrum ---
        nexttile;
        plot(k_rho_vals, avgS_k, 'LineWidth', 2);
        xlabel('k_\rho [\mum^{-1}]','Interpreter','tex','FontSize',axisFontSize);
        ylabel('Magnitude (a.u.)','Interpreter','tex','FontSize',axisFontSize);
        title('Average S(k_\rho)','FontSize',titleFontSize,'FontWeight','bold');
        set(gca,'FontSize',axisFontSize,'YScale','log','XLim',[min(k_rho_vals), max(k_rho_vals)]);
        grid on;

        % --- 3. Angular Spectrum ---
        nexttile;
        plot(theta_vals/pi, avgS_theta, 'LineWidth', 2);
        xlabel('\theta/\pi [rad]','Interpreter','tex','FontSize',axisFontSize);
        ylabel('Magnitude (a.u.)','Interpreter','tex','FontSize',axisFontSize);
        title('Average S(\theta)','FontSize',titleFontSize,'FontWeight','bold');
        set(gca,'FontSize',axisFontSize,'YScale','log','YLim',[1e4, 1e7]);
        ax = gca;
        ax.XMinorGrid = 'on';
        ax.YMinorGrid = 'on';
        grid on;

        drawnow;

        % --- Save figure ---
        saveFigure(fig, ...
            'SaveFileName', opts.SaveFileName, ...
            'SaveDirectory', opts.SaveDirectory, ...
            'SkipSaveFigures', opts.SkipSaveFigures);
    end
end
