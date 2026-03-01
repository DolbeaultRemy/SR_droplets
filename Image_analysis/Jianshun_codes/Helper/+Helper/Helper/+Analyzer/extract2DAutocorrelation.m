function results = extract2DAutocorrelation(od_imgs, scan_parameter_values, options)
%% extract2DAutocorrelation
% Author:       Karthik
% Date:         2025-09-14
% Version:      1.0
%
% Description:
%   Computes 2D autocorrelation g² for a set of real-space images.
%   Returns all g2 maps grouped by unique scan parameter, mean, error
%
% Inputs:
%   od_imgs               - 1xN cell array of 2D image matrices (i × j)
%   scan_parameter_values - array or cell array of scan parameters (similar format as in extractAutocorrelation)
%   options               - options struct containing maximum pixel shift to compute g2 along x and y
%
% Outputs (struct):
%   results.g2_curves             - cell array of 3D autocorrelation matrices [2*max_shift+1 × 2*max_shift+1 × N_reps] per unique scan parameter
%   results.g2_mean               - mean autocorrelation per unique scan parameter
%   results.g2_error              - SEM per unique scan parameter
%   results.scan_parameter_values - unique scan parameter values
%
% Notes:
%   Optional notes, references.

    if isfield(options, 'maximumShift') && ~isempty(options.maximumShift)
        max_shift    = options.maximumShift;
    else
        max_shift    = 5; % [µm]
    end

    pixel_size       = options.pixel_size;
    magnification    = options.magnification;

    % ===== Convert images to 3D array =====
    N_shots = numel(od_imgs);
    [M, N] = size(od_imgs{1});
    img_stack = zeros(M, N, N_shots);
    
    for k = 1:N_shots
        img_stack(:, :, k) = od_imgs{k};
    end
    
    % ===== Determine unique scan parameter values =====
    if isnumeric(scan_parameter_values) && isvector(scan_parameter_values)
        [unique_scan_parameter_values, ~, idx] = unique(scan_parameter_values(:), 'stable');
    elseif iscell(scan_parameter_values)
        params = cell2mat(scan_parameter_values);
        [unique_scan_parameter_values, ~, idx] = unique(params, 'rows', 'stable');
    else
        error('Unsupported format for scan_parameter_values.');
    end
    
    N_unique = size(unique_scan_parameter_values, 1);
    
    % ===== Preallocate outputs =====
    g2_curves = cell(1, N_unique);  % each cell: [2*max_shift+1 × 2*max_shift+1 × N_reps]
    
    % ===== Compute 2D autocorrelation for each unique scan parameter =====
    for i = 1:N_unique
        group_idx  = find(idx == i);            % indices of repetitions for this unique value
        group_data = img_stack(:, :, group_idx);
        N_reps     = length(group_idx);
        
        g2_stack = zeros(2*max_shift+1, 2*max_shift+1, N_reps);
        
        for j = 1:N_reps
            img = group_data(:, :, j);
            [g2_matrix, ~, ~] = Calculator.compute2DAutocorrelation(img, max_shift, pixel_size, magnification);
            g2_stack(:, :, j) = g2_matrix;
        end
        
        % Store all repetitions for this unique scan parameter
        g2_curves{i} = g2_stack;
    end
    
    % ===== Compute mean and SEM per unique value =====
    g2_mean  = cellfun(@(G) mean(G, 3, 'omitnan'), g2_curves, 'UniformOutput', false);
    g2_error = cellfun(@(G) std(G, 0, 3, 'omitnan') ./ sqrt(size(G,3)), g2_curves, 'UniformOutput', false);
    
    % ===== Package results =====
    results = struct();
    results.g2_curves             = g2_curves;                  % raw [2*max_shift+1 × 2*max_shift+1 × N_reps] per unique group
    results.g2_mean               = g2_mean;                    % mean per unique group
    results.g2_error              = g2_error;                   % SEM per unique group
    results.scan_parameter_values = unique_scan_parameter_values;
    results.max_shift             = max_shift;
end
