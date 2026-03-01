function results = conductSpectralAnalysis(od_imgs, scan_parameter_values, options)
%% conductSpectralAnalysis
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Performs Fourier analysis on a set of optical density (OD) images.
%   Computes radial and angular spectral distributions, optionally plots
%   results, and saves figures.
%
% Inputs:
%   od_imgs                  - cell array of OD images
%   scan_parameter_values    - array of scan parameter values corresponding to each image
%   OPTIONS - 
%   saveDirectory            - base directory to save results
%   skipSaveFigures          - skip saving plots
%   skipPreprocessing        - skip preprocessing of images before FFT
%   skipMasking              - skip masking of OD images
%   skipIntensityThresholding- skip thresholding of intensity
%   skipBinarization         - skip binarization of OD images
%   skipNormalization        - skip normalization when plotting angular spectrum
%   skipLivePlot             - skip live plotting of figures
%   pixel_size               - physical pixel size of camera sensor (m)
%   magnification            - imaging magnification
%   zoom_size                - number of pixels to crop around FFT center
%   k_min, k_max             - min/max wavenumber for spectral contrast
%   N_angular_bins           - number of angular bins for S(θ)
%   Angular_Threshold        - threshold parameter for angular spectrum
%   Angular_Sigma            - Gaussian smoothing width for angular spectrum
%   theta_min, theta_max     - angular range for radial spectrum integration
%   N_radial_bins            - number of radial bins for S(k)
%   radial_window_size        - window size for smoothing radial spectrum
%   font                     - font name for plots
%
% Outputs:
%   results                  - struct containing spectra and analysis results
%   Figures (if enabled) are saved into:
%       [saveDirectory]/Results/SpectralAnalysisSavedFigures/
%
% Notes:
%   Optional notes, references.
    
    %% ===== Unpack struct arguments =====
    pixel_size                = options.pixel_size;
    magnification             = options.magnification;
    zoom_size                 = options.zoom_size;
    k_min                     = options.k_min;
    k_max                     = options.k_max;
    N_angular_bins            = options.N_angular_bins;
    Angular_Threshold         = options.Angular_Threshold;
    Angular_Sigma             = options.Angular_Sigma;
    theta_min                 = options.theta_min;
    theta_max                 = options.theta_max;
    N_radial_bins             = options.N_radial_bins;
    radial_window_size        = options.Radial_WindowSize;
    skipNormalization         = options.skipNormalization;
    skipPreprocessing         = options.skipPreprocessing;
    skipMasking               = options.skipMasking;
    skipIntensityThresholding = options.skipIntensityThresholding;
    skipBinarization          = options.skipBinarization;
    skipLivePlot              = options.skipLivePlot;
    skipSaveFigures           = options.skipSaveFigures;
    saveDirectory             = options.saveDirectory;
    font                      = options.font;

    %% ===== Initialization =====
    N_shots                  = length(od_imgs);            % total number of images
    fft_imgs                 = cell(1, N_shots);           % FFT of each image
    angular_spectral_distribution    = cell(1, N_shots);   % S(θ) angular spectrum
    radial_spectral_contrast = zeros(1, N_shots);          % radial contrast metric
    angular_spectral_weight  = zeros(1, N_shots);          % integrated angular weight
    
    S_theta_all              = cell(1, N_shots);
    S_k_all                  = cell(1, N_shots);
    S_k_smoothed_all         = cell(1, N_shots);
    S_theta_norm_all         = cell(1, N_shots);
    PS_all                   = cell(1, N_shots);           % 2D FFT power spectrum |F(kx,ky)|^2

    % Prepare folder to save figures
    if ~skipSaveFigures
        saveFolder = fullfile(saveDirectory, 'Results', 'SavedFigures', 'SpectralAnalysis');
        if ~exist(saveFolder, 'dir')
            mkdir(saveFolder);
        end
    end

    % --- Handle units: allow single string or cell array ---
    if ischar(options.scanParameterUnits) || isstring(options.scanParameterUnits)
        unitList = {char(options.scanParameterUnits)};   % wrap single unit in cell
    else
        unitList = options.scanParameterUnits;          % assume cell array
    end
    
    %% ===== Main loop over images =====
    for k = 1:N_shots
        IMG         = od_imgs{k};
        
        % Skip FFT if image is empty or has low intensity
        if ~(max(IMG(:)) > 1)
            IMGFFT  = NaN(size(IMG));
        else
            % Compute FFT with optional preprocessing
            [IMGFFT, ~] = Calculator.computeFourierTransform(IMG, skipPreprocessing, skipMasking, skipIntensityThresholding, skipBinarization);
        end
    
        % Image size
        [Ny, Nx]    = size(IMG);
    
        % Real-space pixel size (meters)
        dx          = pixel_size / magnification;
        dy          = dx;  % assume square pixels
    
        % Real-space axes in µm
        x           = ((1:Nx) - ceil(Nx/2)) * dx * 1E6;
        y           = ((1:Ny) - ceil(Ny/2)) * dy * 1E6;
    
        % Reciprocal space increments
        dvx         = 1 / (Nx * dx);
        dvy         = 1 / (Ny * dy);
    
        % Frequency axes
        vx          = (-floor(Nx/2):ceil(Nx/2)-1) * dvx;
        vy          = (-floor(Ny/2):ceil(Ny/2)-1) * dvy;
    
        % Wavenumber axes (µm⁻¹)
        kx_full     = 2 * pi * vx * 1E-6;
        ky_full     = 2 * pi * vy * 1E-6;
    
        % Crop FFT image around center
        mid_x       = floor(Nx/2);
        mid_y       = floor(Ny/2);
        fft_imgs{k} = IMGFFT(mid_y-zoom_size:mid_y+zoom_size, mid_x-zoom_size:mid_x+zoom_size);
    
        % Crop wavenumber axes to match cropped FFT
        kx          = kx_full(mid_x - zoom_size : mid_x + zoom_size);
        ky          = ky_full(mid_y - zoom_size : mid_y + zoom_size);
    
        %% ===== Spectral analysis =====
        % Angular spectrum
        [theta_vals, S_theta]       = Calculator.computeAngularSpectralDistribution(fft_imgs{k}, kx, ky, k_min, k_max, N_angular_bins, Angular_Threshold, Angular_Sigma, []);
    
        % Radial spectrum
        [k_rho_vals, S_k]           = Calculator.computeRadialSpectralDistribution(fft_imgs{k}, kx, ky, theta_min, theta_max, N_radial_bins);
    
        % Smooth radial spectrum
        S_k_smoothed                = movmean(S_k, radial_window_size);
    
        % Store results
        angular_spectral_distribution{k}    = S_theta;
        radial_spectral_contrast(k)         = Calculator.computeRadialSpectralContrast(k_rho_vals, S_k_smoothed, k_min, k_max);
    
        % Normalize angular spectrum and compute weight
        S_theta_norm                = S_theta / max(S_theta);
        angular_spectral_weight(k)  = trapz(theta_vals, S_theta_norm);
    
        % Store results
        S_theta_all{k}              = S_theta;
        S_k_all{k}                  = S_k;
        S_k_smoothed_all{k}         = S_k_smoothed;
        S_theta_norm_all{k}         = S_theta_norm;
        PS_all{k}                   = abs(fft_imgs{k}).^2;

        %% ===== Plotting =====
        if ~skipLivePlot
            figure(1); clf
            set(gcf,'Position',[500 100 1000 800])
            tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
            % OD image
            ax1 = nexttile;
            imagesc(x, y, IMG)
            hold on;
            Helper.drawODOverlays(x(1), y(1), x(end), y(end));
            Helper.drawODOverlays(x(end), y(1), x(1), y(end));
            hold off;
            axis equal tight;
            set(gca, 'FontSize', 14, 'YDir', 'normal')
            colormap(ax1, Colormaps.inferno());
            hcb = colorbar;
            ylabel(hcb, 'Optical Density', 'Rotation', -90, 'FontSize', 14, 'FontName', font);
            xlabel('x (\mum)', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            ylabel('y (\mum)', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            title('OD Image', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'tex', 'FontName', font);
    
            % Annotate scan parameter
            % Extract parameter row for this shot
            if iscell(scan_parameter_values)
                % Multi-parameter scan stored as cell array of row vectors
                param_row = scan_parameter_values{k};
            else
                % Numeric vector / matrix
                param_row = scan_parameter_values(k,:);
            end
            
            % --- Ensure units list is long enough ---
            if numel(unitList) < numel(param_row)
                unitList(end+1:numel(param_row)) = {''};  % pad with empty units
            end
            
            % --- Place one text object per parameter ---
            xPos = 0.975;  % normalized x-position (right-aligned)
            yPos = 0.975;  % starting y-position (top)
            yStep = 0.075;  % vertical spacing between multiple parameters
            
            for j = 1:numel(param_row)
                [unitSuffix, txtInterpreter] = getUnitInfo(unitList{j});
                text(xPos, yPos - (j-1)*yStep, sprintf('%.2f%s', param_row(j), unitSuffix), ...
                    'Color', 'white', 'FontWeight', 'bold', 'FontSize', 14, ...
                    'Interpreter', txtInterpreter, 'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            end
            
            % FFT power spectrum
            ax2 = nexttile;
            imagesc(kx, ky, log(1 + PS_all{k}));
            hold on;
            Helper.drawPSOverlays(kx, ky, k_min, k_max)
            
            % Restrict axes strictly to image limits
            xlim([min(kx), max(kx)]);
            ylim([min(ky), max(ky)]);
            axis image; % preserves aspect ratio
            
            set(gca, 'FontSize', 14, 'YDir', 'normal')
            xlabel('k_x [\mum^{-1}]', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            ylabel('k_y [\mum^{-1}]', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            title('Power Spectrum - S(k_x,k_y)', 'Interpreter', 'tex', ...
                  'FontSize', 16, 'FontWeight', 'bold', 'FontName', font);
            colorbar;
            colormap(ax2, Colormaps.coolwarm());
            
            
            % Radial distribution
            nexttile;
            plot(k_rho_vals, S_k_smoothed, 'LineWidth', 2);
            set(gca, 'FontSize', 14, 'YScale', 'log', 'XLim', [min(k_rho_vals), max(k_rho_vals)]);
            xlabel('k_\rho [\mum^{-1}]', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            ylabel('Magnitude (a.u.)', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            title('Radial Spectral Distribution - S(k_\rho)', 'Interpreter', 'tex', ...
                'FontSize', 16, 'FontWeight', 'bold', 'FontName', font);
            grid on;

            % Angular distribution
            nexttile;
            if ~skipNormalization
                plot(theta_vals/pi, S_theta_norm, 'LineWidth', 2);
                set(gca, 'FontSize', 14, 'YLim', [0, 1]);
            else
                plot(theta_vals/pi, S_theta, 'LineWidth', 2);
                set(gca, 'FontSize', 14, 'YScale', 'log', 'YLim', [1E4, 1E7]);
            end
            xlabel('\theta/\pi [rad]', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            ylabel('Magnitude (a.u.)', 'Interpreter', 'tex', 'FontSize', 14, 'FontName', font);
            title('Angular Spectral Distribution - S(\theta)', 'Interpreter', 'tex', ...
                'FontSize', 16, 'FontWeight', 'bold', 'FontName', font);
            grid on;
            ax = gca;
            ax.MinorGridLineStyle = ':';
            ax.MinorGridColor = [0.7 0.7 0.7];
            ax.MinorGridAlpha = 0.5;
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
        end

        %% ===== Save figures =====
        if ~skipSaveFigures
            fileNamePNG = fullfile(saveFolder, sprintf('fft_analysis_img_%03d.png', k));
            print(gcf, fileNamePNG, '-dpng', '-r100');
        elseif ~skipLivePlot
            pause(0.5);
        end
    end
    
     % Package results into struct
    results = struct();
    results.kx                                = kx;
    results.ky                                = ky;
    results.PS_all                            = PS_all;
    results.theta_vals                        = theta_vals;
    results.S_theta_all                       = S_theta_all;
    results.k_rho_vals                        = k_rho_vals;
    results.S_k_all                           = S_k_all;
    results.angular_spectral_distribution     = angular_spectral_distribution;
    results.S_k_smoothed_all                  = S_k_smoothed_all;
    results.radial_spectral_contrast          = radial_spectral_contrast;
    results.S_theta_norm_all                  = S_theta_norm_all;
    results.angular_spectral_weight           = angular_spectral_weight;
    
end

%% === Local helper function ===
function [unitSuffix, txtInterpreter] = getUnitInfo(u)
    switch lower(u)
        case {'degrees','deg','°'}
            unitSuffix     = '^\circ';
            txtInterpreter = 'tex';
        case {'gauss','g'}
            unitSuffix     = ' G';
            txtInterpreter = 'none';
        otherwise
            unitSuffix     = '';
            txtInterpreter = 'none';
    end
end