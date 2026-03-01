classdef TwoGaussian2DModel < handle
%% TwoGaussian2DModel
% Author: Jianshun Gao
% Date: 2025-09-12
% Version: 1.0
%
% Description:
%   This class provides methods to model, fit, and analyze 2D datasets using
%   a sum of two Gaussian distributions. It supports:
%     - Initial parameter estimation from 2D data.
%     - Nonlinear least-squares fitting to extract Gaussian parameters.
%     - Access and modification of parameter values and bounds.
%
% Properties:
%   - params: Structure storing parameter values and bounds.
%   - fitResult: MATLAB fit object after fitting.
%   - gof: Goodness-of-fit structure from the fit.
%
% Methods:
%   - guess(data, x, y): Estimate initial parameters from 2D data.
%   - fit(data, x, y): Fit the two-Gaussian model to data.
%   - getParamValue(paramName): Retrieve the current value of a parameter.
%   - setParamValue(paramName, value): Update the value of a parameter.
%   - setParamBounds(paramName, minVal, maxVal): Set parameter bounds.
%
% Usage:
%   obj = TwoGaussian2DModel();                       % Create object
%   params = obj.guess(data, x, y);                  % Estimate parameters
%   [fitResult, gof] = obj.fit(data, x, y);          % Perform fitting
%   ampA = obj.getParamValue('A_amplitude');          % Access parameter
%   obj.setParamValue('B_sigmax', 5.0);               % Modify parameter
%   obj.setParamBounds('A_sigmay', 1.0, 10.0);       % Set bounds
%
% Notes:
%   - All parameter names must match those defined in the params structure.
%   - This class is intended for 2D datasets where two Gaussian peaks are present.

    properties
        params;     % Structure storing parameter values, bounds, etc.
        fitResult;  % MATLAB fit object after fitting
        gof;        % Goodness-of-fit metrics
    end
    
    methods
        function obj = TwoGaussian2DModel()
            % Constructor (empty)
        end
        
        function params = guess(obj, data, x, y, varargin)
            % Estimate initial parameters for two Gaussian peaks
            % Usage: obj.guess(data, x, y)
            % Optional parameter:
            %   'negative' - boolean, allow negative amplitude (default: false)
            
            p = inputParser;
            addParameter(p, 'negative', false);
            parse(p, varargin{:});
            
            % Simple peak estimation
            [maxVal, maxIdx] = max(data(:));
            [row, col] = ind2sub(size(data), maxIdx);
            
            % Create parameter structure with values, min, and max
            params = struct();
            
            % Amplitude parameters
            params.A_amplitude.value = maxVal;
            params.A_amplitude.min = 0;
            params.A_amplitude.max = 2 * maxVal;
            
            params.B_amplitude.value = maxVal / 2;
            params.B_amplitude.min = 0;
            params.B_amplitude.max = 2 * maxVal;
            
            % Center positions
            params.A_centerx.value = x(col);
            params.A_centerx.min = min(x);
            params.A_centerx.max = max(x);
            
            params.A_centery.value = y(row);
            params.A_centery.min = min(y);
            params.A_centery.max = max(y);
            
            params.B_centerx.value = x(round(end/2));
            params.B_centerx.min = min(x);
            params.B_centerx.max = max(x);
            
            params.B_centery.value = y(round(end/2));
            params.B_centery.min = min(y);
            params.B_centery.max = max(y);
            
            % Standard deviations
            sigmax_range = max(x) - min(x);
            sigmay_range = max(y) - min(y);
            
            params.A_sigmax.value = sigmax_range / 10;
            params.A_sigmax.min = sigmax_range / 100;
            params.A_sigmax.max = sigmax_range / 2;
            
            params.A_sigmay.value = sigmay_range / 10;
            params.A_sigmay.min = sigmay_range / 100;
            params.A_sigmay.max = sigmay_range / 2;
            
            params.B_sigmax.value = sigmax_range / 5;
            params.B_sigmax.min = sigmax_range / 100;
            params.B_sigmax.max = sigmax_range / 2;
            
            params.B_sigmay.value = sigmay_range / 5;
            params.B_sigmay.min = sigmay_range / 100;
            params.B_sigmay.max = sigmay_range / 2;
            
            obj.params = params;
        end
        
        function [fitResult, gof] = fit(obj, data, x, y, varargin)
            % Fit the data to a sum of two 2D Gaussian distributions
            % Usage: [fitResult, gof] = obj.fit(data, x, y)
            % Optional parameter:
            %   'params' - struct with initial guesses and bounds
            
            p = inputParser;
            addParameter(p, 'params', []);
            parse(p, varargin{:});
            
            if isempty(p.Results.params) && isempty(obj.params)
                obj.guess(data, x, y, varargin{:});
            elseif ~isempty(p.Results.params)
                obj.params = p.Results.params;
            end
            
            % Prepare fitting data
            [X, Y] = meshgrid(x, y);
            xData = X(:);
            yData = Y(:);
            zData = data(:);
            
            % Define two-Gaussian function
            twoGauss = @(A_amplitude, B_amplitude, A_centerx, A_centery, ...
                B_centerx, B_centery, A_sigmax, A_sigmay, B_sigmax, B_sigmay, x, y) ...
                A_amplitude * exp(-((x-A_centerx).^2/(2*A_sigmax^2) + (y-A_centery).^2/(2*A_sigmay^2))) + ...
                B_amplitude * exp(-((x-B_centerx).^2/(2*B_sigmax^2) + (y-B_centery).^2/(2*B_sigmay^2)));
            
            % Create MATLAB fittype
            ft = fittype(twoGauss, 'independent', {'x','y'}, 'dependent', 'z');
            
            % Fit options
            options = fitoptions(ft);
            
            % Parameter order
            paramOrder = {'A_amplitude', 'B_amplitude', 'A_centerx', 'A_centery', ...
                          'B_centerx', 'B_centery', 'A_sigmax', 'A_sigmay', 'B_sigmax', 'B_sigmay'};
            
            % Build StartPoint, Lower, Upper
            startPoint = zeros(1, length(paramOrder));
            lowerBounds = zeros(1, length(paramOrder));
            upperBounds = zeros(1, length(paramOrder));
            for i = 1:length(paramOrder)
                paramName = paramOrder{i};
                startPoint(i) = obj.params.(paramName).value;
                lowerBounds(i) = obj.params.(paramName).min;
                upperBounds(i) = obj.params.(paramName).max;
            end
            
            options.StartPoint = startPoint;
            options.Lower = lowerBounds;
            options.Upper = upperBounds;
            
            % Perform the fit
            [obj.fitResult, obj.gof] = fit([xData, yData], zData, ft, options);
            fitResult = obj.fitResult;
            gof = obj.gof;
        end
        
        function paramValue = getParamValue(obj, paramName)
            % Get the value of a specified parameter
            if isfield(obj.params, paramName)
                paramValue = obj.params.(paramName).value;
            else
                error('Parameter %s does not exist', paramName);
            end
        end
        
        function setParamValue(obj, paramName, value)
            % Set the value of a specified parameter
            if isfield(obj.params, paramName)
                obj.params.(paramName).value = value;
            else
                error('Parameter %s does not exist', paramName);
            end
        end
        
        function setParamBounds(obj, paramName, minVal, maxVal)
            % Set the bounds of a specified parameter
            if isfield(obj.params, paramName)
                obj.params.(paramName).min = minVal;
                obj.params.(paramName).max = maxVal;
            else
                error('Parameter %s does not exist', paramName);
            end
        end
    end

end
