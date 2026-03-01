function results = extractCustomCorrelation(angular_spectral_distribution, scan_parameter_values, N_shots, N_angular_bins)
%% extractCustomCorrelation
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Extracts correlation of a single (highest) peak with possible secondary peak (50-70°)
%   Computes correlations and cumulants grouped by unique scan parameter.
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
    [unique_scan_parameter_values, ~, idx] = unique(scan_parameter_values(:), 'stable');
elseif iscell(scan_parameter_values)
    params = cell2mat(scan_parameter_values);
    [unique_scan_parameter_values, ~, idx] = unique(params, 'rows', 'stable');
else
    error('Unsupported format for scan_parameter_values.');
end

N_unique = size(unique_scan_parameter_values, 1);

% ===== Angular settings =====
angle_range     = 180;
angle_per_bin   = angle_range / N_angular_bins;
max_peak_bin    = round(180 / angle_per_bin);
window_size     = 10;
angle_threshold = 100;

% ===== Preallocate result arrays based on unique values =====
mean_max_g2_values                        = zeros(1, N_unique);
skew_max_g2_angle_values                  = zeros(1, N_unique);
var_max_g2_values                         = zeros(1, N_unique);
fourth_order_cumulant_max_g2_angle_values = zeros(1, N_unique);
max_g2_all_per_scan_parameter_value       = cell(1, N_unique);
std_error_g2_values                       = zeros(1, N_unique);

% ===== Compute correlations and cumulants per unique parameter =====
for i = 1:N_unique
    group_idx  = find(idx == i);          % all repetitions for this unique value
    group_data = delta_nkr_all(group_idx, :);
    N_reps     = size(group_data, 1);

    g2_values = zeros(1, N_reps);

    for j = 1:N_reps
        profile = group_data(j, :);

        % Find highest peak in 0–180° range
        restricted_profile = profile(1:max_peak_bin);
        [~, peak_idx_rel] = max(restricted_profile);
        peak_idx          = peak_idx_rel;
        peak_angle        = (peak_idx - 1) * angle_per_bin;

        % Determine offsets for secondary peak correlation
        if peak_angle < angle_threshold
            offsets = round(50 / angle_per_bin) : round(70 / angle_per_bin);
        else
            offsets = -round(70 / angle_per_bin) : -round(50 / angle_per_bin);
        end

        ref_window = mod((peak_idx - window_size):(peak_idx + window_size) - 1, N_angular_bins) + 1;
        ref        = profile(ref_window);

        correlations = zeros(size(offsets));
        for k_off = 1:length(offsets)
            shifted_idx = mod(peak_idx + offsets(k_off) - 1, N_angular_bins) + 1;
            sec_window  = mod((shifted_idx - window_size):(shifted_idx + window_size) - 1, N_angular_bins) + 1;
            sec         = profile(sec_window);
            correlations(k_off) = mean(ref .* sec) / mean(ref.^2);
        end

        [max_corr, ~] = max(correlations);
        g2_values(j)  = max_corr;
    end

    % Store all repetitions for this unique parameter
    max_g2_all_per_scan_parameter_value{i} = g2_values;

    % Compute cumulants
    kappa = Calculator.computeCumulants(g2_values(:), 4);

    mean_max_g2_values(i)                        = kappa(1);
    var_max_g2_values(i)                         = kappa(2);
    skew_max_g2_angle_values(i)                  = kappa(3);
    fourth_order_cumulant_max_g2_angle_values(i) = kappa(4);

    N_eff = sum(~isnan(g2_values));
    std_error_g2_values(i) = sqrt(kappa(2)) / sqrt(N_eff);
end

% ===== Package results into struct =====
results = struct();
results.mean_max_g2                         = mean_max_g2_values;
results.var_max_g2                          = var_max_g2_values;
results.skew_max_g2_angle                   = skew_max_g2_angle_values;
results.fourth_order_cumulant_max_g2        = fourth_order_cumulant_max_g2_angle_values;
results.std_error_g2                        = std_error_g2_values;
results.max_g2_all_per_scan_parameter_value = max_g2_all_per_scan_parameter_value;
results.scan_parameter_values               = unique_scan_parameter_values;
end