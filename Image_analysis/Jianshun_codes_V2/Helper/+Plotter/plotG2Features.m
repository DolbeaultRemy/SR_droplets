function plotG2Features(results, varargin)
%% plotG2Features
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Plot Fourier amplitudes, contrast, and symmetry fractions.
%
% Notes:
%   Optional notes, references.

    % --- Parse name-value pairs ---
    p = inputParser;
    addParameter(p, 'Title', 'Analysis Results', @(x) ischar(x) || isstring(x));
    addParameter(p, 'XLabel', 'Scan Parameter', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName', 'Arial', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'FigNum', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'SkipSaveFigures', false, @islogical);
    addParameter(p, 'SaveFileName', 'AllFeatures.fig', @ischar);
    addParameter(p, 'SaveDirectory', pwd, @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    T = results.summaries;
    x = T.scanParamVal;

    % --- Extract Fourier amplitudes, contrast, symmetry fractions ---
    A2m = T.mean_A2; A2e = T.Fun_A2;
    A4m = T.mean_A4; A4e = T.Fun_A4;
    A6m = T.mean_A6; A6e = T.Fun_A6;
    Nrm = T.mean_Contrast; Nrme = T.Fun_Contrast;
    S2m = T.mean_S2; S2e = T.Fun_S2;
    Q4m = T.mean_Q4; Q4e = T.Fun_Q4;
    H6m = T.mean_H6; H6e = T.Fun_H6;

    % --- Color palette: same color for each symmetry type ---
    colors = [...
        0.8500 0.3250 0.0980;  % Stripe (red/orange) -> A2/S2
        0.0000 0.4470 0.7410;  % Square (blue) -> A4/Q4
        0.4660 0.6740 0.1880]; % Hexatic (green) -> A6/H6

    % --- Create figure ---
    if isempty(opts.FigNum), fig = figure; else, fig = figure(opts.FigNum); end
    clf(fig);
    set(fig,'Color','w','Position',[100 100 950 750]);
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

    % --- Top-left tile: Fourier amplitudes ---
    ax1 = nexttile; hold on;
    errorbar(x, A2m, A2e,'-o','LineWidth',1.5,'Color',colors(1,:),'MarkerFaceColor',colors(1,:));
    errorbar(x, A4m, A4e,'-s','LineWidth',1.5,'Color',colors(2,:),'MarkerFaceColor',colors(2,:));
    errorbar(x, A6m, A6e,'-^','LineWidth',1.5,'Color',colors(3,:),'MarkerFaceColor',colors(3,:));
    xlabel(ax1,opts.XLabel,'FontName',opts.FontName,'FontSize',opts.FontSize);
    ylabel(ax1,'Amplitude','FontName',opts.FontName,'FontSize',opts.FontSize);
    title(ax1,'Harmonic Amplitudes','FontName',opts.FontName,'FontSize',opts.FontSize+2);
    legend(ax1,'|A2|','|A4|','|A6|','Location','best'); grid(ax1,'on');
    set(ax1,'FontName',opts.FontName,'FontSize',opts.FontSize);
    
    % --- Top-right tile: Symmetry fractions ---
    ax2 = nexttile; hold on;
    errorbar(x, S2m, S2e,'-o','LineWidth',1.5,'Color',colors(1,:),'MarkerFaceColor',colors(1,:));
    errorbar(x, Q4m, Q4e,'-s','LineWidth',1.5,'Color',colors(2,:),'MarkerFaceColor',colors(2,:));
    errorbar(x, H6m, H6e,'-^','LineWidth',1.5,'Color',colors(3,:),'MarkerFaceColor',colors(3,:));
    xlabel(ax2,opts.XLabel,'FontName',opts.FontName,'FontSize',opts.FontSize);
    ylabel(ax2,'Symmetry Fraction','FontName',opts.FontName,'FontSize',opts.FontSize);
    ylim(ax2,[0 1]); grid(ax2,'on');
    title(ax2,'Symmetry Fractions','FontName',opts.FontName,'FontSize',opts.FontSize+2);
    legend(ax2,'S2','Q4','H6','Location','best');
    set(ax2,'FontName',opts.FontName,'FontSize',opts.FontSize);
    
    % --- Bottom-left tile: Contrast (mean of individual curves + mean curve) ---
    ax3 = nexttile; 
    % Compute mean ± SEM of individual curve contrasts
    N_params = numel(results.features_group);
    meanContrast_ind = zeros(1,N_params);
    semContrast_ind  = zeros(1,N_params);
    for i = 1:N_params
        FT = results.features_group{i};
        contrasts = FT.Contrast;
        meanContrast_ind(i) = mean(contrasts,'omitnan');
        semContrast_ind(i)  = std(contrasts,0,1,'omitnan')/sqrt(numel(contrasts));
    end
    errorbar(x, meanContrast_ind, semContrast_ind,'-d','LineWidth',1.5,'Color',[0.2 0.6 0.2],'MarkerFaceColor',[0.2 0.6 0.2]);
    hold on;
    % Plot mean-curve contrast on top
    errorbar(x, Nrm, Nrme,'-d','LineWidth',1.5,'Color',[0.6 0.2 0.8],'MarkerFaceColor',[0.6 0.2 0.8]);
    xlabel(ax3,opts.XLabel,'FontName',opts.FontName,'FontSize',opts.FontSize);
    ylabel(ax3,'Contrast','FontName',opts.FontName,'FontSize',opts.FontSize);
    title(ax3,'g^2 Curve Contrast','FontName',opts.FontName,'FontSize',opts.FontSize+2);
    grid(ax3,'on'); set(ax3,'FontName',opts.FontName,'FontSize',opts.FontSize);
    legend(ax3,'Individual Curves Mean ± SEM','Mean Curve','Location','southwest');

    % --- Bottom-right tile: Max-peak location vs scan parameter ---
    ax4 = nexttile; hold on;
    % Extract max peak for each scan parameter and convert to degrees
    maxPeaks = rad2deg(cellfun(@(s) s.MaxPeakAngle, results.meanCurve));
    plot(x, maxPeaks,'-o','LineWidth',1.5,'Color',[0.8 0.3 0.0],'MarkerFaceColor',[0.8 0.3 0.0]);
    ylim(ax4, [10, 90])
    xlabel(ax4,opts.XLabel,'FontName',opts.FontName,'FontSize',opts.FontSize);
    ylabel(ax4,'\theta (degrees)','Interpreter','tex','FontSize',opts.FontSize);
    title(ax4,'Max peak in mean g^2 curve','FontName',opts.FontName,'FontSize',opts.FontSize+2);
    grid(ax4,'on'); set(ax4,'FontName',opts.FontName,'FontSize',opts.FontSize);

    %% --- Save figure ---
    if ~opts.SkipSaveFigures
        if ~exist(opts.SaveDirectory,'dir'), mkdir(opts.SaveDirectory); end
        savefig(fig, fullfile(opts.SaveDirectory, opts.SaveFileName));
    end
end