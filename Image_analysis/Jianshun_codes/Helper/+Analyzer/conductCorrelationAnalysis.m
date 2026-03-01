function results = conductCorrelationAnalysis(od_imgs, scan_parameter_values, options)
%% conductCorrelationAnalysis
% Author:       Karthik
% Date:         2025-09-14
% Version:      1.0
%
% Description:
%   Computes 2D autocorrelation g²(Δx,Δy) on OD images.
%   Extracts radial and angular distributions using Calculator functions.
%   Optionally plots results and saves figures.
%
% Inputs:
%   od_imgs                  - cell array of OD images
%   scan_parameter_values    - array of scan parameter values
%   options - struct with fields:
%       saveDirectory        - base directory to save results
%       skipSaveFigures      - skip saving plots
%       skipLivePlot         - skip live plotting
%       pixel_size           - physical pixel size of camera sensor (m)
%       magnification        - imaging magnification
%       maximumShift         - maximum pixel shift for g²
%       font                 - font name for plots
%
% Outputs:
%   results - struct containing g² maps, radial and angular distributions
%
% Notes:
%   Optional notes, references.

    %% ===== Unpack struct arguments =====
    pixel_size       = options.pixel_size;
    magnification    = options.magnification;
    skipLivePlot     = options.skipLivePlot;
    skipSaveFigures  = options.skipSaveFigures;
    saveDirectory    = options.saveDirectory;
    font             = options.font;

    radial_theta       = options.Radial_Theta;
    radial_window_size = options.Radial_WindowSize;
    N_angular_bins     = options.N_angular_bins;
    r_min              = options.Radial_Minimum;
    r_max              = options.Radial_Maximum;

    if isfield(options, 'maximumShift') && ~isempty(options.maximumShift)
        maximumShift = options.maximumShift;
    else
        maximumShift = 5; % [µm]
    end

    % --- Handle units ---
    if ischar(options.scanParameterUnits) || isstring(options.scanParameterUnits)
        unitList = {char(options.scanParameterUnits)};
    else
        unitList = options.scanParameterUnits;
    end

    %% ===== Initialization =====
    N_shots     = numel(od_imgs);
    g2_matrices = cell(1, N_shots);
    g2_radial   = cell(1, N_shots);
    g2_angular  = cell(1, N_shots);
    r_vals      = cell(1, N_shots);
    theta_vals  = cell(1, N_shots);

    dx          = pixel_size / magnification;
    shifts      = -maximumShift:maximumShift;

    if ~skipSaveFigures
        saveFolder = fullfile(saveDirectory, 'Results', 'SavedFigures', 'AutocorrAnalysis');
        if ~exist(saveFolder, 'dir')
            mkdir(saveFolder);
        end
    end

    %% ===== Loop over images =====
    for k = 1:N_shots
        IMG = od_imgs{k};

        % Compute g² in physical units
        [g2_matrix, dx_phys, dy_phys] = Calculator.compute2DAutocorrelation( ...
            IMG, maximumShift, pixel_size, magnification);
        g2_matrices{k} = g2_matrix;
        
        % Extract radial profile (within angular window)
        [r_vals{k}, g2_radial{k}] = Calculator.computeRadialCorrelation( ...
            g2_matrix, dx_phys, dy_phys, radial_theta);

        % Extract angular profile (within radial band)
        [theta_vals{k}, g2_angular{k}] = Calculator.computeAngularCorrelation( ...
            g2_matrix, dx_phys, dy_phys, r_min, r_max, N_angular_bins);
        
        % Smooth radial profile
        g2_radial_smoothed = movmean(g2_radial{k}, radial_window_size);

        
        %% ===== Plotting =====
        if ~skipLivePlot
            figure(1); clf
            set(gcf,'Position',[500 100 1000 800])
            tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

            % OD image
            ax1 = nexttile;
            imagesc(IMG);
            axis equal tight;
            set(gca,'FontSize',14,'YDir','normal');
            colormap(ax1, Colormaps.inferno());
            hcb = colorbar;
            ylabel(hcb,'Optical Density','Rotation',-90,'FontSize',14,'FontName',font);
            xlabel('x [px]','FontSize',14,'FontName',font);
            ylabel('y [px]','FontSize',14,'FontName',font);
            title('OD Image','FontSize',16,'FontWeight','bold','FontName',font);

            % Annotate scan parameter
            if iscell(scan_parameter_values)
                param_row = scan_parameter_values{k};
            else
                param_row = scan_parameter_values(k,:);
            end
            if numel(unitList) < numel(param_row)
                unitList(end+1:numel(param_row)) = {''};
            end
            xPos = 0.975; yPos = 0.975; yStep = 0.075;
            for j = 1:numel(param_row)
                [unitSuffix, txtInterpreter] = getUnitInfo(unitList{j});
                text(xPos, yPos-(j-1)*yStep, sprintf('%.2f%s',param_row(j),unitSuffix), ...
                    'Color','white','FontWeight','bold','FontSize',14, ...
                    'Interpreter',txtInterpreter,'Units','normalized', ...
                    'HorizontalAlignment','right','VerticalAlignment','top');
            end

            % g² map
            ax2 = nexttile;
            imagesc(shifts, shifts, g2_matrix);
            axis equal tight;
            set(gca,'FontSize',14,'YDir','normal');
            colormap(ax2, Colormaps.coolwarm());
            colorbar;
            xlabel('\Deltax (\mum)','FontSize',14,'FontName',font);
            ylabel('\Deltay (\mum)','FontSize',14,'FontName',font);
            title('Autocorrelation g_2(\Deltax,\Deltay)','FontSize',16,'FontWeight','bold','FontName',font);

            % Radial distribution
            ax3 = nexttile;
            plot(r_vals{k}, g2_radial_smoothed, 'LineWidth',2);
            set(gca,'FontSize',14);
            set(gca, 'FontSize', 14, 'XLim', [min(r_vals{k}), max(r_vals{k})], 'YLim', [0, 1]);
            xlabel('r [\mum]','Interpreter','tex','FontSize',14,'FontName',font);
            ylabel('g_2(r)','FontSize',14,'FontName',font);
            title('Radial Distribution','FontSize',16,'FontWeight','bold','FontName',font);
            grid on;

            % Angular distribution
            ax4 = nexttile;
            plot(theta_vals{k}/pi, g2_angular{k}, 'LineWidth',2);
            set(gca,'FontSize', 14, 'YLim', [0, 1]);
            xlabel('\theta/\pi [rad]','Interpreter','tex','FontSize',14,'FontName',font);
            ylabel('g_2(\theta)','FontSize',14,'FontName',font);
            title('Angular Distribution','FontSize',16,'FontWeight','bold','FontName',font);
            grid on;
            ax = gca;
            ax.MinorGridLineStyle = ':';
            ax.MinorGridColor = [0.7 0.7 0.7];
            ax.MinorGridAlpha = 0.5;
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
        end

        %% ===== Save figures =====
        if ~skipSaveFigures && ~skipLivePlot
            fileNamePNG = fullfile(saveFolder, sprintf('g2_analysis_img_%03d.png', k));
            print(gcf, fileNamePNG, '-dpng','-r100');
        elseif ~skipLivePlot
            pause(0.5);
        end
    end

    %% ===== Package results =====
    results = struct();
    results.g2_matrices           = g2_matrices;
    results.g2_radial             = g2_radial;
    results.g2_angular            = g2_angular;
    results.r_vals                = r_vals;
    results.theta_vals            = theta_vals;
    results.scan_parameter_values = unique(scan_parameter_values);

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
