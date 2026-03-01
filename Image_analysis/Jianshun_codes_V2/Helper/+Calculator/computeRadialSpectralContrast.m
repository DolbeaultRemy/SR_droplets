function contrast = computeRadialSpectralContrast(k_rho_vals, S_k_smoothed, k_min, k_max)
%% computeRadialSpectralContrast
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Computes the ratio of the peak in S_k_smoothed within [k_min, k_max] to the value at (or near) k = 0.
%
% Notes:
%   Optional notes, references.

    % Ensure inputs are column vectors
    k_rho_vals     = k_rho_vals(:);
    S_k_smoothed   = S_k_smoothed(:);

    % Step 1: Find index of k ≈ 0
    [~, idx_k0] = min(abs(k_rho_vals));  % Closest to zero
    S_k0 = S_k_smoothed(idx_k0);

    % Step 2: Find indices in specified k-range
    in_range = (k_rho_vals >= k_min) & (k_rho_vals <= k_max);
    
    if ~any(in_range)
        warning('No values found in the specified k-range. Returning NaN.');
        contrast = NaN;
        return;
    end

    % Step 3: Find peak value in the specified k-range
    S_k_peak = max(S_k_smoothed(in_range));

    % Step 4: Compute contrast
    contrast = S_k_peak / S_k0;

end