function [k_rho_vals, S_radial] = computeRadialSpectralDistribution(IMGFFT, kx, ky, thetamin, thetamax, num_bins)
%% computeRadialSpectralDistribution
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Computes the ratio of the peak in S_k_smoothed within [k_min, k_max] to the value at (or near) k = 0.
%
% Inputs:
%   IMGFFT    : 2D FFT image (fftshifted and cropped)
%   kx, ky    : 1D physical wavenumber axes [μm⁻¹] matching FFT size
%   thetamin  : Minimum angle (in radians)
%   thetamax  : Maximum angle (in radians)
%   num_bins  : Number of radial bins
%
% Outputs:
%   k_rho_vals:
%   S_radial  :  
% 
% Notes:
%   Optional notes, references.

    [KX, KY] = meshgrid(kx, ky);
    K_rho = sqrt(KX.^2 + KY.^2);
    Theta = atan2(KY, KX);

    if thetamin < thetamax
        angle_mask = (Theta >= thetamin) & (Theta <= thetamax);
    else
        angle_mask = (Theta >= thetamin) | (Theta <= thetamax);
    end

    power_spectrum = abs(IMGFFT).^2;

    r_min = min(K_rho(angle_mask));
    r_max = max(K_rho(angle_mask));
    r_edges = linspace(r_min, r_max, num_bins + 1);
    k_rho_vals = 0.5 * (r_edges(1:end-1) + r_edges(2:end));
    S_radial = zeros(1, num_bins);

    for i = 1:num_bins
        r_low = r_edges(i);
        r_high = r_edges(i + 1);
        radial_mask = (K_rho >= r_low) & (K_rho < r_high);
        full_mask = radial_mask & angle_mask;
        S_radial(i) = sum(power_spectrum(full_mask));
    end
end