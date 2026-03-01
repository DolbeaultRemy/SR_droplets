function results = extractAutocorrelation(theta_values, angular_spectral_distribution, scan_parameter_values, N_shots, N_angular_bins)
%% extractAutocorrelation
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Computes angular autocorrelation g² for a set of angular spectral distribution.
%   Returns all g2 curves grouped by unique scan parameter, mean, error
%
% Notes:
%   Optional notes, references.

% ===== Convert spectral distributions to matrix =====
delta_nkr_all = zeros(N_shots, N_angular_bins);
for k = 1:N_shots
    delta_nkr_all(k, :) = angular_spectral_distribution{k};
end

% ===== Determine unique scan parameter values =====
if isnumeric(scan_parameter_values) && isvector(scan_parameter_values)
    % Single parameter case
    [unique_scan_parameter_values, ~, idx] = unique(scan_parameter_values(:), 'stable');
elseif iscell(scan_parameter_values)
    % Multi-parameter case (each row is a parameter set)
    params = cell2mat(scan_parameter_values);
    [unique_scan_parameter_values, ~, idx] = unique(params, 'rows', 'stable');
else
    error('Unsupported format for scan_parameter_values.');
end

N_unique = size(unique_scan_parameter_values, 1);

% ===== Preallocate outputs =====
g2_curves = cell(1, N_unique);   % each cell: [N_reps × N_angular_bins]

% ===== Compute g²(θ) for each unique scan parameter =====
for i = 1:N_unique
    group_idx  = find(idx == i);            % indices of repetitions for this unique value
    group_data = delta_nkr_all(group_idx, :);
    N_reps     = length(group_idx);
    
    g2_matrix = zeros(N_reps, N_angular_bins);
    
    for j = 1:N_reps
        profile = group_data(j, :);
        for dtheta = 0:N_angular_bins-1
            profile_shifted = circshift(profile, -dtheta, 2);
            g2_matrix(j, dtheta+1) = mean(profile .* profile_shifted) / mean(profile.^2);
        end
    end
    
    % Store all repetitions for this unique scan parameter
    g2_curves{i} = g2_matrix;
end

% ===== Compute mean and SEM per unique value =====
g2_mean  = cellfun(@(G) mean(G, 1, 'omitnan'), g2_curves, 'UniformOutput', false);
g2_error = cellfun(@(G) std(G, 0, 1, 'omitnan') ./ sqrt(size(G,1)), g2_curves, 'UniformOutput', false);

% ===== Package results =====
results = struct();
results.g2_curves             = g2_curves;                  % raw [N_reps × N_angular_bins] per unique group
results.g2_mean               = g2_mean;                    % mean per unique group
results.g2_error              = g2_error;                   % SEM per unique group
results.theta_values          = theta_values;
results.scan_parameter_values = unique_scan_parameter_values;
end
