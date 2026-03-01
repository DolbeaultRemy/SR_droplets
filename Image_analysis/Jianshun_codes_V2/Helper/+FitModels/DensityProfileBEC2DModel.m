classdef DensityProfileBEC2DModel < handle
%% DensityProfileBEC2DModel
% Author: Jianshun Gao
% Date: 2025-09-12
% Version: 1.0
%
% Description:
%   This class provides methods to model, fit, and analyze 2D Bose-Einstein 
%   condensate (BEC) density profiles, including both condensed (Thomas-Fermi)
%   and thermal (polylog) components.
%
% Properties:
%   - params: Structure storing parameter values and bounds.
%   - fitResult: MATLAB fit object after fitting.
%   - gof: Goodness-of-fit structure from the fit.
%
% Methods:
%   - cal_cond_frac(X, Y): Compute condensate fraction from fitted profile.
%   - return_atom_number(X, Y): Compute total, BEC, and thermal atom numbers.
%   - return_temperature(tof, omega): Estimate temperature from thermal cloud width.
%
% Static Methods:
%   - ThomasFermi_2d: 2D Thomas-Fermi parabolic profile.
%   - polylog: Series approximation of polylogarithm function.
%   - polylog2_2d: 2D thermal cloud profile using polylog(2) function.
%   - density_profile_BEC_2d: Full 2D density (BEC + thermal), optionally rotated.
%   - density_1d: 1D slice of the density profile combining BEC and thermal parts.
%
% Usage:
%   obj       = DensityProfileBEC2DModel();
%   atom_n    = obj.return_atom_number(X, Y);        % Atom numbers
%   cond_frac = obj.cal_cond_frac(X, Y);             % Condensate fraction
%   T         = obj.return_temperature(tof, omega);  % Temperature
%
% Notes:
%   - All static methods can be called independently without creating an object.

    properties
        % Conversion factors and default constants
        fwhm_factor   = 2*sqrt(2*log(2));  % FWHM factor for Gaussian
        height_factor = 1/(2*pi);          % Height normalization factor
        prefix        = '';                % Optional parameter prefix
        tiny          = 1e-15;             % Small number to avoid division by zero
        params;                            % Struct for parameter definitions
        fitResult;                         % Result of 2D fit
        cond_frac;                         % Condensate fraction
        gof;                               % Goodness-of-fit structure
        atom_n_conv   = 144;               % Conversion factor for atom number
        pre_check     = false;             % Flag for pre-fit check
        post_check    = false;             % Flag for post-fit check
        is_debug      = false;             % Flag for debug plots
    end
    
    methods
        function obj = DensityProfileBEC2DModel(prefix)
            % Constructor
            if nargin > 0
                obj.prefix = prefix;
            end
            obj.set_paramhints_prefix();
        end
        
        function set_paramhints_prefix(obj)
            % Initialize parameter hints and expressions
            obj.params = struct();
            
            % Minimum values
            obj.params.([obj.prefix 'amp_bec']).min = 0;
            obj.params.([obj.prefix 'amp_th']).min = 0;
            obj.params.([obj.prefix 'x0_bec']).min = 0;
            obj.params.([obj.prefix 'y0_bec']).min = 0;
            obj.params.([obj.prefix 'x0_th']).min = 0;
            obj.params.([obj.prefix 'y0_th']).min = 0;
            obj.params.([obj.prefix 'sigmax_bec']).min = 0;
            obj.params.([obj.prefix 'sigmay_bec']).min = 0;
            obj.params.([obj.prefix 'sigma_th']).min = 0;
            
            obj.params.([obj.prefix 'rot_angle']).min = -90;
            obj.params.([obj.prefix 'rot_angle']).max = 90;
            
            % Expressions for derived parameters
            obj.params.([obj.prefix 'atom_number_bec']).expr = ...
                [obj.prefix 'amp_bec / 5 * 2 * pi * ' obj.prefix 'sigmax_bec * ' obj.prefix 'sigmay_bec'];
            obj.params.([obj.prefix 'atom_number_th']).expr = ...
                [obj.prefix 'amp_th * 2 * pi * 1.20206 / 1.643 * ' obj.prefix 'sigma_th * ' obj.prefix 'sigma_th'];
            obj.params.([obj.prefix 'condensate_fraction']).expr = ...
                [obj.prefix 'atom_number_bec / (' obj.prefix 'atom_number_bec + ' obj.prefix 'atom_number_th)'];
        end
        
        function params = guess(obj, data, x, y, varargin)
            % Estimate initial parameters for BEC + thermal cloud fit.
            % Performs 1D slicing along the shorter axis of the image and initial amplitude guesses.
            
            p = inputParser;
            addParameter(p, 'negative', false);
            addParameter(p, 'pureBECThreshold', 0.5);
            addParameter(p, 'noBECThThreshold', 0.0);
            addParameter(p, 'rot_angle', 0);
            addParameter(p, 'vary_rot', false);
            addParameter(p, 'is_debug', false);
            addParameter(p, 'pre_check', false);
            addParameter(p, 'post_check', false);
            parse(p, varargin{:});
            
            obj.is_debug = p.Results.is_debug;
            obj.pre_check = p.Results.pre_check;
            obj.post_check = p.Results.post_check;
            
            x_width = length(unique(x));
            y_width = length(unique(y));
            x_1d = linspace(x(1), x(end), x_width);
            y_1d = linspace(y(1), y(end), y_width);
            
            data_2d = reshape(data, y_width, x_width)';
            
            % Rotate image if needed
            rot_angle = p.Results.rot_angle;
            if rot_angle ~= 0
                data_2d = imrotate(data_2d, rot_angle, 'bilinear', 'crop');
            end
            
            % Debug plot of input data
            if obj.is_debug
                figure;
                imagesc(x_1d, y_1d, data_2d');
                axis equal tight;
                colorbar;
                colormap('jet');
                title('Input Data');
                xlabel('x-axis');
                ylabel('y-axis');
            end
            
            % Binarize and locate BEC
            thresh = obj.calc_thresh(data_2d, 0.5);
            center_pix = obj.calc_cen_pix(thresh);
            center = obj.center_pix_conv(center_pix, x_1d, y_1d);
            BEC_width_guess = obj.guess_BEC_width(thresh, center_pix);
        
            if obj.is_debug
                figure;
                imagesc(x_1d, y_1d, thresh');
                hold on;
                plot(center(1), center(2), 'gx', 'MarkerSize', 25, 'LineWidth', 2);
                axis equal tight;
                colorbar;
                colormap('jet');
                title(sprintf('Binarized Image (BEC Width: x=%.0f, y=%.0f pixels)', BEC_width_guess(1), BEC_width_guess(2)));
                xlabel('x-axis');
                ylabel('y-axis');
            end
            
            % 1D slicing along the shorter axis
            if BEC_width_guess(1) < BEC_width_guess(2)
                if obj.is_debug
                    disp('x-axis is shorter, performing 1D fit along x-axis');
                end
                s_width_ind = 1;
                x_fit = x_1d;
                slice_range = round(center_pix(2) - BEC_width_guess(2)/2) : round(center_pix(2) + BEC_width_guess(2)/2);
                X_guess = sum(data_2d(:, slice_range), 2) / length(slice_range);
            else
                if obj.is_debug
                    disp('y-axis is shorter, performing 1D fit along y-axis');
                end
                s_width_ind = 2;
                x_fit = y_1d;
                slice_range = round(center_pix(1) - BEC_width_guess(1)/2) : round(center_pix(1) + BEC_width_guess(1)/2);
                X_guess = sum(data_2d(slice_range, :), 1) / length(slice_range);
            end
            
            max_val = max(X_guess);
                    
            % Construct initial parameter struct
            params_1d = struct();
            params_1d.x0_bec.value = center(s_width_ind);
            params_1d.x0_bec.min = center(s_width_ind) - 10;
            params_1d.x0_bec.max = center(s_width_ind) + 10;
            
            params_1d.x0_th.value = center(s_width_ind);
            params_1d.x0_th.min = center(s_width_ind) - 10;
            params_1d.x0_th.max = center(s_width_ind) + 10;
            
            params_1d.amp_bec.value = 0.5 * max_val;
            params_1d.amp_bec.min = 0;
            params_1d.amp_bec.max = 1.3 * max_val;
            
            params_1d.amp_th.value = 0.5 * max_val;
            params_1d.amp_th.min = 0;
            params_1d.amp_th.max = 1.3 * max_val;
            
            % params_1d.deltax.value = 3 * BEC_width_guess(s_width_ind);
            % params_1d.deltax.min = 0;
            % params_1d.deltax.max = max(x_width, y_width);
            
            params_1d.sigma_bec.value = BEC_width_guess(s_width_ind) / 1.22;
            params_1d.sigma_bec.min = 0;
            params_1d.sigma_bec.max = 2 * BEC_width_guess(s_width_ind);
            
            % params_1d.sigma_th.value = 3 * BEC_width_guess(1);
            params_1d.sigma_th.value = 0.632 * params_1d.sigma_bec.value + 0.518 * 3 * BEC_width_guess(s_width_ind);
            params_1d.sigma_th.min = 0;
            params_1d.sigma_th.max = 3 * (0.632 * params_1d.sigma_bec.value + 0.518 * 3 * BEC_width_guess(s_width_ind));
            % params_1d.sigma_th.expr = '0.632*sigma_bec + 0.518*deltax';
        
            % Perform 1D bimodal fit
            [fitResult_1d, gof_1d] = fit_1d_bimodal(x_fit, X_guess, params_1d);
            
            % Extract fit coefficients
            bval_1d = coeffvalues(fitResult_1d);
            paramNames_1d = coeffnames(fitResult_1d);
            
            for i = 1:length(paramNames_1d)
                bval_1d_struct.(paramNames_1d{i}) = bval_1d(i);
            end
            
            if obj.is_debug
                disp('Initial 1D fit parameters:');
                disp(bval_1d);
                figure;
                plot(x_fit, X_guess, 'b-', 'LineWidth', 2);
                hold on;
                plot(x_fit, obj.density_1d(x_fit, bval_1d.x0_bec, bval_1d.x0_th, ...
                    bval_1d.amp_bec, bval_1d.amp_th, bval_1d.sigma_bec, bval_1d.sigma_th), 'r-', 'LineWidth', 2);
                legend('Data', 'Fit');
                if s_width_ind == 1
                    xlabel('x-axis'); title('1D Fit along x-axis');
                else
                    xlabel('y-axis'); title('1D Fit along y-axis');
                end
                ylabel('Intensity');
            end
            
            % Scale amplitudes
            blurred_data = imgaussfilt(data_2d, 1);
            amp_conv_1d_2d = max(blurred_data(:)) / (bval_1d_struct.amp_bec + bval_1d_struct.amp_th);
            max_val = max(data_2d(:));
            
            % Create parameter struct
            params = struct();
            
            % Pre-check: decide if image is pure BEC or pure thermal cloud based on 1D fit result
            if bval_1d_struct.amp_th / bval_1d_struct.amp_bec > 7 && obj.pre_check
                if obj.is_debug
                    disp('Image seems to be pure thermal cloud (based on 1D fit amplitudes)');
                end
                
                params.([obj.prefix 'amp_bec']).value = 0;
                params.([obj.prefix 'amp_bec']).vary = false;
                
                params.([obj.prefix 'amp_th']).value = amp_conv_1d_2d * bval_1d_struct.amp_th;
                params.([obj.prefix 'amp_th']).max = 1.3 * max_val;
                params.([obj.prefix 'amp_th']).vary = true;
                
                params.([obj.prefix 'x0_bec']).value = 1;
                params.([obj.prefix 'x0_bec']).vary = false;
                
                params.([obj.prefix 'y0_bec']).value = 1;
                params.([obj.prefix 'y0_bec']).vary = false;
                
                params.([obj.prefix 'x0_th']).value = center(1);
                params.([obj.prefix 'x0_th']).min = center(1) - 10;
                params.([obj.prefix 'x0_th']).max = center(1) + 10;
                params.([obj.prefix 'x0_th']).vary = true;
                
                params.([obj.prefix 'y0_th']).value = center(2);
                params.([obj.prefix 'y0_th']).min = center(2) - 10;
                params.([obj.prefix 'y0_th']).max = center(2) + 10;
                params.([obj.prefix 'y0_th']).vary = true;
                
                params.([obj.prefix 'sigmax_bec']).value = 1;
                params.([obj.prefix 'sigmax_bec']).vary = false;
                
                params.([obj.prefix 'sigmay_bec']).value = 1;
                params.([obj.prefix 'sigmay_bec']).vary = false;
                
                params.([obj.prefix 'sigma_th']).value = bval_1d_struct.sigma_th;
                params.([obj.prefix 'sigma_th']).min = 0;
                params.([obj.prefix 'sigma_th']).max = max(x_width, y_width);
                params.([obj.prefix 'sigma_th']).vary = true;
                
            elseif bval_1d_struct.amp_bec / bval_1d_struct.amp_th > 10 && obj.pre_check
                if obj.is_debug
                    disp('Image seems to be pure BEC (based on 1D fit amplitudes)');
                end
                
                params.([obj.prefix 'amp_bec']).value = amp_conv_1d_2d * bval_1d_struct.amp_bec;
                params.([obj.prefix 'amp_bec']).max = 1.3 * max_val;
                params.([obj.prefix 'amp_bec']).vary = true;
                
                params.([obj.prefix 'amp_th']).value = 0;
                params.([obj.prefix 'amp_th']).vary = false;
                
                params.([obj.prefix 'x0_bec']).value = center(1);
                params.([obj.prefix 'x0_bec']).min = center(1) - 10;
                params.([obj.prefix 'x0_bec']).max = center(1) + 10;
                params.([obj.prefix 'x0_bec']).vary = true;
                
                params.([obj.prefix 'y0_bec']).value = center(2);
                params.([obj.prefix 'y0_bec']).min = center(2) - 10;
                params.([obj.prefix 'y0_bec']).max = center(2) + 10;
                params.([obj.prefix 'y0_bec']).vary = true;
                
                params.([obj.prefix 'x0_th']).value = 1;
                params.([obj.prefix 'x0_th']).vary = false;
                
                params.([obj.prefix 'y0_th']).value = 1;
                params.([obj.prefix 'y0_th']).vary = false;
                
                params.([obj.prefix 'sigma_th']).value = 1;
                params.([obj.prefix 'sigma_th']).vary = false;
                
                if s_width_ind == 1
                    params.([obj.prefix 'sigmax_bec']).value = bval_1d_struct.sigma_bec;
                    params.([obj.prefix 'sigmax_bec']).max = 2 * BEC_width_guess(1);
                    params.([obj.prefix 'sigmax_bec']).vary = true;
                    
                    params.([obj.prefix 'sigmay_bec']).value = BEC_width_guess(2) / 1.22;
                    params.([obj.prefix 'sigmay_bec']).max = 2 * BEC_width_guess(2);
                    params.([obj.prefix 'sigmay_bec']).vary = true;
                else
                    params.([obj.prefix 'sigmax_bec']).value = BEC_width_guess(1) / 1.22;
                    params.([obj.prefix 'sigmax_bec']).max = 2 * BEC_width_guess(1);
                    params.([obj.prefix 'sigmax_bec']).vary = true;
                    
                    params.([obj.prefix 'sigmay_bec']).value = bval_1d_struct.sigma_bec;
                    params.([obj.prefix 'sigmay_bec']).max = 2 * BEC_width_guess(2);
                    params.([obj.prefix 'sigmay_bec']).vary = true;
                end
                
            else
                % Normal bimodal fit parameters
                params.([obj.prefix 'amp_bec']).value = amp_conv_1d_2d * bval_1d_struct.amp_bec;
                params.([obj.prefix 'amp_bec']).min = 0;
                params.([obj.prefix 'amp_bec']).max = 1.3 * max_val;
                params.([obj.prefix 'amp_bec']).vary = true;
                
                params.([obj.prefix 'amp_th']).value = amp_conv_1d_2d * bval_1d_struct.amp_th;
                params.([obj.prefix 'amp_th']).min = 0;
                params.([obj.prefix 'amp_th']).max = 1.3 * max_val;
                params.([obj.prefix 'amp_th']).vary = true;
                
                params.([obj.prefix 'x0_bec']).value = center(1);
                params.([obj.prefix 'x0_bec']).min = center(1) - 10;
                params.([obj.prefix 'x0_bec']).max = center(1) + 10;
                params.([obj.prefix 'x0_bec']).vary = true;
                
                params.([obj.prefix 'y0_bec']).value = center(2);
                params.([obj.prefix 'y0_bec']).min = center(2) - 10;
                params.([obj.prefix 'y0_bec']).max = center(2) + 10;
                params.([obj.prefix 'y0_bec']).vary = true;
                
                params.([obj.prefix 'x0_th']).value = center(1);
                params.([obj.prefix 'x0_th']).min = center(1) - 10;
                params.([obj.prefix 'x0_th']).max = center(1) + 10;
                params.([obj.prefix 'x0_th']).vary = true;
                
                params.([obj.prefix 'y0_th']).value = center(2);
                params.([obj.prefix 'y0_th']).min = center(2) - 10;
                params.([obj.prefix 'y0_th']).max = center(2) + 10;
                params.([obj.prefix 'y0_th']).vary = true;
                
                params.([obj.prefix 'sigma_th']).value = bval_1d_struct.sigma_th;
                params.([obj.prefix 'sigma_th']).min = 0;
                params.([obj.prefix 'sigma_th']).max = max(x_width, y_width);
                params.([obj.prefix 'sigma_th']).vary = true;
                
                if s_width_ind == 1
                    params.([obj.prefix 'sigmax_bec']).value = bval_1d_struct.sigma_bec;
                    params.([obj.prefix 'sigmax_bec']).min = 0;
                    params.([obj.prefix 'sigmax_bec']).max = 2 * BEC_width_guess(1);
                    params.([obj.prefix 'sigmax_bec']).vary = true;
                    
                    params.([obj.prefix 'sigmay_bec']).value = BEC_width_guess(2) / 1.22;
                    params.([obj.prefix 'sigmay_bec']).min = 0;
                    params.([obj.prefix 'sigmay_bec']).max = 2 * BEC_width_guess(2);
                    params.([obj.prefix 'sigmay_bec']).vary = true;
                else
                    params.([obj.prefix 'sigmax_bec']).value = BEC_width_guess(1) / 1.22;
                    params.([obj.prefix 'sigmax_bec']).min = 0;
                    params.([obj.prefix 'sigmax_bec']).max = 2 * BEC_width_guess(1);
                    params.([obj.prefix 'sigmax_bec']).vary = true;
                    
                    params.([obj.prefix 'sigmay_bec']).value = bval_1d_struct.sigma_bec;
                    params.([obj.prefix 'sigmay_bec']).min = 0;
                    params.([obj.prefix 'sigmay_bec']).max = 2 * BEC_width_guess(2);
                    params.([obj.prefix 'sigmay_bec']).vary = true;
                end
            end
            
            params.([obj.prefix 'rot_angle']).value = rot_angle;
            params.([obj.prefix 'rot_angle']).min = rot_angle - 30;
            params.([obj.prefix 'rot_angle']).max = rot_angle + 30;
            params.([obj.prefix 'rot_angle']).vary = p.Results.vary_rot;
            
            if obj.is_debug
                disp('Initial parameters:');
                disp(params);
            end
            
            obj.params = params;
        end
        
        function [fitResult, gof] = fit(obj, data, x, y, varargin)
            data = double(data);
            
            % Perform fitting
            if isempty(obj.params)
                obj.guess(data, x, y, varargin{:});
            end
        
            % Prepare fitting data
            [X, Y] = meshgrid(x, y);
            xyData = [X(:), Y(:)];
            zData = double(data(:));
        
            % Create fit options
            options = fitoptions('Method', 'NonlinearLeastSquares');
        
            if obj.params.([obj.prefix 'rot_angle']).vary
                % Define parameter order
                paramOrder = {[obj.prefix 'amp_bec'], [obj.prefix 'amp_th'], ...
                             [obj.prefix 'x0_bec'], [obj.prefix 'y0_bec'], ...
                             [obj.prefix 'x0_th'], [obj.prefix 'y0_th'], ...
                             [obj.prefix 'sigmax_bec'], [obj.prefix 'sigmay_bec'], ...
                             [obj.prefix 'sigma_th'], [obj.prefix 'rot_angle']};
        
                % Create StartPoint, Lower, and Upper vectors
                startPoint = zeros(1, length(paramOrder));
                lowerBounds = zeros(1, length(paramOrder));
                upperBounds = zeros(1, length(paramOrder));
        
                for i = 1:length(paramOrder)
                    paramName = paramOrder{i};
                    startPoint(i) = obj.params.(paramName).value;
                    lowerBounds(i) = obj.params.(paramName).min;
                    upperBounds(i) = obj.params.(paramName).max;
                end
        
                % Set fitting options
                options.StartPoint = startPoint;
                options.Lower = lowerBounds;
                options.Upper = upperBounds;
        
                % Define fit type
                ft = fittype(@(amp_bec, amp_th, x0_bec, y0_bec, x0_th, y0_th, ...
                    sigmax_bec, sigmay_bec, sigma_th, rot_angle, x, y) ...
                    obj.density_profile_BEC_2d(x, y, amp_bec, amp_th, x0_bec, y0_bec, ...
                    x0_th, y0_th, sigmax_bec, sigmay_bec, sigma_th, rot_angle), ...
                    'independent', {'x', 'y'}, 'dependent', 'z');
            else
                % Define parameter order
                paramOrder = {[obj.prefix 'amp_bec'], [obj.prefix 'amp_th'], ...
                             [obj.prefix 'x0_bec'], [obj.prefix 'y0_bec'], ...
                             [obj.prefix 'x0_th'], [obj.prefix 'y0_th'], ...
                             [obj.prefix 'sigmax_bec'], [obj.prefix 'sigmay_bec'], ...
                             [obj.prefix 'sigma_th']};
        
                % Create StartPoint, Lower, and Upper vectors
                startPoint = zeros(1, length(paramOrder));
                lowerBounds = zeros(1, length(paramOrder));
                upperBounds = zeros(1, length(paramOrder));
        
                for i = 1:length(paramOrder)
                    paramName = paramOrder{i};
                    startPoint(i) = obj.params.(paramName).value;
                    lowerBounds(i) = obj.params.(paramName).min;
                    upperBounds(i) = obj.params.(paramName).max;
                end
        
                % Set fitting options
                options.StartPoint = startPoint;
                options.Lower = lowerBounds;
                options.Upper = upperBounds;
        
                % Define fit type
                ft = fittype(@(amp_bec, amp_th, x0_bec, y0_bec, x0_th, y0_th, ...
                    sigmax_bec, sigmay_bec, sigma_th, x, y) ...
                    obj.density_profile_BEC_2d(x, y, amp_bec, amp_th, x0_bec, y0_bec, ...
                    x0_th, y0_th, sigmax_bec, sigmay_bec, sigma_th, 0), ...
                    'independent', {'x', 'y'}, 'dependent', 'z');
            end
        
            % Perform fitting
            [obj.fitResult, obj.gof] = fit(xyData, zData, ft, options);
            fitResult = obj.fitResult;
            gof = obj.gof;
        
            % Post-processing check
            if obj.post_check
                bval = coeffvalues(obj.fitResult);
                paramNames = coeffnames(obj.fitResult);
        
                % Extract parameter values
                for i = 1:length(paramNames)
                    eval([paramNames{i} ' = bval(i);']);
                end
        
                % Calculate number of atoms around the BEC
                tf_fit = obj.ThomasFermi_2d(xyData(:,1), xyData(:,2), x0_bec, y0_bec, amp_bec, sigmax_bec, sigmay_bec);
                tf_fit_2 = obj.ThomasFermi_2d(xyData(:,1), xyData(:,2), x0_bec, y0_bec, amp_bec, 1.5*sigmax_bec, 1.5*sigmay_bec);
        
                mask = tf_fit > 0;
                mask_2 = tf_fit_2 > 0;
        
                N_c = sum(zData(mask & ~mask_2));
                N_a = obj.atom_n_conv * N_c;
        
                % If too few atoms are found around the BEC, refit (BEC only)
                if N_a < 6615
                    if obj.is_debug
                        disp('No thermal cloud detected, performing BEC-only fit');
                    end
        
                    % Update parameters
                    obj.params.([obj.prefix 'amp_th']).value = 0;
                    obj.params.([obj.prefix 'amp_th']).vary = false;
        
                    obj.params.([obj.prefix 'x0_th']).value = 1;
                    obj.params.([obj.prefix 'x0_th']).vary = false;
        
                    obj.params.([obj.prefix 'y0_th']).value = 1;
                    obj.params.([obj.prefix 'y0_th']).vary = false;
        
                    obj.params.([obj.prefix 'sigma_th']).value = 1;
                    obj.params.([obj.prefix 'sigma_th']).vary = false;
        
                    % Refit
                    [obj.fitResult, obj.gof] = fit(xyData, zData, ft, options);
                    fitResult = obj.fitResult;
                    gof = obj.gof;
                end
            end
        
            % Calculate condensate fraction
            obj.cond_frac = obj.cal_cond_frac(X, Y);
        end
        
        function thresh = calc_thresh(obj, data, thresh_val, sigma)
            % Binarize image
            if nargin < 3
                thresh_val = 0.3;
            end
            if nargin < 4
                sigma = 0.4;
            end
        
            blurred = imgaussfilt(data, sigma);
            thresh = blurred < max(blurred(:)) * thresh_val;
            thresh = double(~thresh); % Invert and convert to double precision
        end
        
        function center_pix = calc_cen_pix(obj, thresh)
            % Calculate the center of the binarized image
            [X, Y] = size(thresh);
        
            thresh = thresh / sum(thresh(:));
        
            % Edge distributions
            dx = sum(thresh, 2);
            dy = sum(thresh, 1);
        
            % Expectation values
            center_pix = [sum(dx .* (1:X)'), sum(dy .* (1:Y))];
        end
        
        function center = center_pix_conv(obj, center_pix, x, y)
            % Convert pixel center to coordinate center
            center = [x(round(center_pix(1))), y(round(center_pix(2)))];
        end
        
        function BEC_width_guess = guess_BEC_width(obj, thresh, center)
            % Guess BEC width
            [X, Y] = size(thresh);
        
            BEC_width_guess = [sum(thresh(:, round(center(2)))), sum(thresh(round(center(1)), :))];
        
            for i = 1:2
                if BEC_width_guess(i) <= 0
                    BEC_width_guess(i) = 1;
                end
            end
        end
        
        function cond_frac = cal_cond_frac(obj, X, Y)
            % Calculate condensate fraction
            bval = coeffvalues(obj.fitResult);
            paramNames = coeffnames(obj.fitResult);
            
            % Extract parameter values
            for i = 1:length(paramNames)
                eval([paramNames{i} ' = bval(i);']);
            end
        
            if ~obj.params.([obj.prefix 'rot_angle']).vary
                rot_angle = 0;
            end
            
            tf_fit = obj.ThomasFermi_2d(X, Y, x0_bec, y0_bec, amp_bec, sigmax_bec, sigmay_bec);
            fit_total = obj.density_profile_BEC_2d(X, Y, amp_bec, amp_th, x0_bec, y0_bec, x0_th, y0_th, sigmax_bec, sigmay_bec, sigma_th, rot_angle);
            
            N_bec = sum(tf_fit(:));
            N_ges = sum(fit_total(:));
            cond_frac = N_bec / N_ges;
        end
        
        function atom_n = return_atom_number(obj, X, Y, is_print)
            % Calculate atom number
            if nargin < 4
                is_print = true;
            end
            
            bval = coeffvalues(obj.fitResult);
            paramNames = coeffnames(obj.fitResult);
            
            % Extract parameter values
            for i = 1:length(paramNames)
                eval([paramNames{i} ' = bval(i);']);
            end
            
            tf_fit = obj.ThomasFermi_2d(X, Y, x0_bec, y0_bec, amp_bec, sigmax_bec, sigmay_bec);
            th_fit = obj.polylog2_2d(X, Y, x0_th, y0_th, amp_th, sigma_th, sigma_th);
            
            N_bec = obj.atom_n_conv * sum(tf_fit(:));
            N_th = obj.atom_n_conv * sum(th_fit(:));
            N = N_bec + N_th;
            frac = N_bec / N;
            
            if is_print
                fprintf('\nAtom numbers:\n');
                fprintf('   N_bec: %.0f\n', N_bec);
                fprintf('   N_th: %.0f\n', N_th);
                fprintf('   N_total: %.0f\n', N);
                fprintf('   Condensate fraction: %.2f %%\n', frac * 100);
            end
            
            atom_n = struct('N', N, 'N_bec', N_bec, 'N_th', N_th, 'cond_f', frac);
        end
        
        function T = return_temperature(obj, tof, omg, is_print, eff_pix)
            % Calculate temperature
            if nargin < 3
                omg = [];
            end
            if nargin < 4
                is_print = true;
            end
            if nargin < 5
                eff_pix = 2.472e-6;
            end
            
            bval = coeffvalues(obj.fitResult);
            paramNames = coeffnames(obj.fitResult);
            
            % Extract parameter values
            for i = 1:length(paramNames)
                eval([paramNames{i} ' = bval(i);']);
            end
            
            R_th = sigma_th * eff_pix * sqrt(2);
            
            % Physical constants
            u = 1.66053906660e-27; % Atomic mass unit
            k = 1.380649e-23;      % Boltzmann constant
            
            if isempty(omg)
                T = R_th^2 * 164 * u / k / tof^2;
            else
                T = R_th^2 * 164 * u / k / (1/omg^2 + tof^2);
            end
            
            if is_print
                fprintf('Temperature: %.2f nK\n', T * 1e9);
            end
        end
    end
    
    methods (Static)

        function res = ThomasFermi_2d(x, y, centerx, centery, amplitude, sigmax, sigmay)
            % Thomas-Fermi distribution function
            % tiny = 1e-15;
            res = (1 - ((x - centerx) / sigmax).^2 - ((y - centery) / sigmay).^2);
            res(res < 0) = 0; 
            res = amplitude * res.^(3/2);
            % res = amplitude * 5/(2*pi) / max(tiny, sigmax * sigmay) .* (res > 0) .* res;
        end
        
        function res = polylog(power, numerator)
            % Polylogarithm function approximation
            order = 20;
            dataShape = size(numerator);
            numerator = repmat(numerator(:), 1, order);
            numerator = numerator .^ repmat(1:order, prod(dataShape), 1);
            
            denominator = repmat((1:order), prod(dataShape), 1);
            data = numerator ./ (denominator .^ power);
            
            res = sum(data, 2);
            res = reshape(res, dataShape);
        end
        
        function res = polylog2_2d(x, y, centerx, centery, amplitude, sigmax, sigmay)
            % 2D polylogarithm function
            % tiny = 1e-15;
            arg = exp(-((x - centerx).^2 / (2 * sigmax^2)) - ((y - centery).^2 / (2 * sigmay^2)));
            res = amplitude / (1.643) .* FitModels.DensityProfileBEC2DModel.polylog(2, arg);
        end
        
        function res = density_profile_BEC_2d(x, y, amp_bec, amp_th, x0_bec, y0_bec, x0_th, y0_th, sigmax_bec, sigmay_bec, sigma_th, rot_angle)
            % BEC density profile function
            if nargin < 12
                rot_angle = 0;
            end
            
            % Rotate coordinates (if needed)
            if rot_angle ~= 0
                rot_angle_rad = -rot_angle * pi/180; % Negative sign means clockwise rotation
                
                % Rotate coordinates
                x_rot = x * cos(rot_angle_rad) + y * sin(rot_angle_rad);
                y_rot = -x * sin(rot_angle_rad) + y * cos(rot_angle_rad);
                
                % Rotate BEC center
                x0_bec_rot = x0_bec * cos(rot_angle_rad) + y0_bec * sin(rot_angle_rad);
                y0_bec_rot = -x0_bec * sin(rot_angle_rad) + y0_bec * cos(rot_angle_rad);
                
                % Rotate thermal center
                x0_th_rot = x0_th * cos(rot_angle_rad) + y0_th * sin(rot_angle_rad);
                y0_th_rot = -x0_th * sin(rot_angle_rad) + y0_th * cos(rot_angle_rad);
                
                x = x_rot;
                y = y_rot;
                x0_bec = x0_bec_rot;
                y0_bec = y0_bec_rot;
                x0_th = x0_th_rot;
                y0_th = y0_th_rot;
            end
            
            % Calculate Thomas-Fermi part
            TF_part = FitModels.DensityProfileBEC2DModel.ThomasFermi_2d(x, y, x0_bec, y0_bec, amp_bec, sigmax_bec, sigmay_bec);
            
            % Calculate polylogarithm part
            poly_part = FitModels.DensityProfileBEC2DModel.polylog2_2d(x, y, x0_th, y0_th, amp_th, sigma_th, sigma_th);
            
            % Total sum
            res = TF_part + poly_part;
        end
        
        function res = density_1d(x, x0_bec, x0_th, amp_bec, amp_th, sigma_bec, sigma_th)
            % 1D density profile (Thomas-Fermi + thermal polylog)
            thermal_part = amp_th / 1.643 * polylog_int(exp(-0.5 * (x - x0_th).^2 / sigma_th^2));
            TF_part = amp_bec * (1 - ((x - x0_bec) / sigma_bec).^2);
            TF_part(TF_part < 0) = 0;
            TF_part = TF_part.^(3/2);
            res = thermal_part + TF_part;
        end
    end

end

% Helper function: polylogarithm interpolation

function res = polylog_int(x)
    % Create interpolation table (simplified version)
    x_int = linspace(0, 1.00001, 1000);
    poly_tab = zeros(size(x_int));
    
    for i = 1:length(x_int)
        poly_tab(i) = sum((x_int(i).^(1:20)) ./ (1:20).^2);
    end
    
    % Linear interpolation
    res = interp1(x_int, poly_tab, x, 'linear', 'extrap');
end

function [fitResult, gof] = fit_1d_bimodal(x, y, initialParams)
    % 1D bimodal fitting function
    % Input:
    %   x - independent variable data
    %   y - dependent variable data
    %   initialParams - structure of initial parameters
    % Output:
    %   fitResult - fitting result
    %   gof - goodness-of-fit statistics
    
    % Define 1D bimodal fitting function
    bimodal1d = @(amp_bec, amp_th, x0_bec, x0_th, sigma_bec, sigma_th, x) ...
        FitModels.DensityProfileBEC2DModel.density_1d(x, x0_bec, x0_th, amp_bec, amp_th, sigma_bec, sigma_th);
    paramNames = {'amp_bec', 'amp_th', 'x0_bec', 'x0_th', 'sigma_bec', 'sigma_th'};
    
    % Create fit type
    ft = fittype(bimodal1d, 'independent', 'x', 'dependent', 'y');
    
    % Set fit options
    options = fitoptions(ft);
    
    % Set initial parameters and bounds
    startPoint = zeros(1, length(paramNames));
    lowerBounds = zeros(1, length(paramNames));
    upperBounds = zeros(1, length(paramNames));
    
    for i = 1:length(paramNames)
        paramName = paramNames{i};
        startPoint(i) = initialParams.(paramName).value;
        lowerBounds(i) = initialParams.(paramName).min;
        upperBounds(i) = initialParams.(paramName).max;
    end
    
    options.StartPoint = startPoint;
    options.Lower = lowerBounds;
    options.Upper = upperBounds;
    
    % Perform fitting
    [fitResult, gof] = fit(x(:), y(:), ft, options);
end