function [results, scan_parameter_values, scan_reference_values] = performAnalysis(options)
%% performAnalysis
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Brief description of the script functionality.
%
% Notes:
%   Optional notes, references.

    arguments
        options.scan_parameter (1,:) char
        options.ignore_scan_parameter {mustBeText} = ''
        options.scan_reference_values (1,:) double
        options.cam (1,1) double
        options.angle (1,1) double
        options.center (1,2) double
        options.span (1,2) double
        options.fraction (1,2) double
        options.ImagingMode (1,:) char
        options.PulseDuration (1,1) double
        options.pixel_size (1,1) double
        options.magnification (1,1) double
        options.zoom_size (1,1) double
        options.N_angular_bins (1,1) double
        options.Angular_Threshold (1,1) double
        options.Angular_Sigma (1,1) double
        options.Angular_WindowSize (1,1) double
        options.theta_min (1,1) double
        options.theta_max (1,1) double
        options.N_radial_bins (1,1) double
        options.Radial_Sigma (1,1) double
        options.Radial_WindowSize (1,1) double
        options.k_min (1,1) double
        options.k_max (1,1) double
        options.skipFringeRemoval (1,1) logical
        options.skipUnshuffling (1,1) logical
        options.flipSortOrder (1,1) logical
        options.skipPreprocessing (1,1) logical
        options.skipMasking (1,1) logical
        options.skipIntensityThresholding (1,1) logical
        options.skipBinarization (1,1) logical
        options.skipNormalization (1,1) logical
        options.skipLivePlot (1,1) logical
        options.skipSaveFigures (1,1) logical
        options.skipSaveData (1,1) logical
        options.skipSaveProcessedOD (1,1) logical
        options.skipFullODImagesFolderUse (1,1) logical
        options.showProgressBar (1,1) logical
        options.measurementName (1,:) char
        options.scanParameterUnits {mustBeText} = ''
        options.selectedPath (1,:) char
        options.folderPath (1,:) char
        options.baseDataFolder (1,:) char
        options.saveDirectory (1,:) char
        options.FullODImagesFolder (1,:) char
        options.titleString (1,:) char
        options.font (1,:) char
        options.SAVE_TO_WORKSPACE (1,1) logical
    end

    % Collect OD images
    [od_imgs, scan_parameter_values, scan_reference_values, ~] = Helper.collectODImages(options);

    % Conduct spectral analysis
    fprintf('\n[INFO] Initiating spectral analysis...\n');
    
    spectral_analysis_results           = Analyzer.conductSpectralAnalysis(od_imgs, scan_parameter_values, options);
    
    N_shots                             = length(od_imgs);
    
    % Extract angular correlations
    full_g2_results                     = Analyzer.extractAutocorrelation(...
                                          spectral_analysis_results.theta_vals, ...
                                          spectral_analysis_results.angular_spectral_distribution, ...
                                          scan_parameter_values, N_shots, options.N_angular_bins);
    
    custom_g_results                    = Analyzer.extractCustomCorrelation(...
                                          spectral_analysis_results.angular_spectral_distribution, ...
                                          scan_parameter_values, N_shots, options.N_angular_bins);
    
    fprintf('\n[INFO] Spectral analysis complete!\n');
    
    % Conduct PCA
    fprintf('\n[INFO] Initiating Principal Component Analysis...\n');
    
    pca_results                         = Analyzer.conductPCA(od_imgs);
    
    fprintf('\n[INFO] Principal Component Analysis complete!\n');

    % Lattice Reconstruction
    
    % Package results into struct
    results = struct();
    results.spectral_analysis_results   = spectral_analysis_results;
    results.full_g2_results             = full_g2_results;
    results.custom_g_results            = custom_g_results;
    results.pca_results                 = pca_results;
end