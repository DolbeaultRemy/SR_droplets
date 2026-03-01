function [theta_vals, S_theta] = computeAngularSpectralDistribution(IMGFFT, kx, ky, k_min, k_max, num_bins, threshold, sigma, windowSize)
%% computeAngularSpectralDistribution
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Brief description of the script functionality.
%
% Notes:
%   Optional notes, references.

    % Apply threshold to isolate strong peaks
    IMGFFT(IMGFFT < threshold) = 0;

    % Create wavenumber meshgrid
    [KX, KY] = meshgrid(kx, ky);
    Kmag     = sqrt(KX.^2 + KY.^2);       % radial wavenumber magnitude
    Theta    = atan2(KY, KX);             % range [-pi, pi]

    % Restrict to radial band in wavenumber space
    radial_mask = (Kmag >= k_min) & (Kmag <= k_max);

    % Initialize angular structure factor
    S_theta    = zeros(1, num_bins);
    theta_vals = linspace(0, pi, num_bins);  % only 0 to pi due to symmetry

    % Loop over angular bins
    for i = 1:num_bins
        angle_start = (i - 1) * pi / num_bins;
        angle_end   = i * pi / num_bins;
        angle_mask  = (Theta >= angle_start) & (Theta < angle_end);
        bin_mask    = radial_mask & angle_mask;
        fft_angle   = IMGFFT .* bin_mask;
        S_theta(i)  = sum(sum(abs(fft_angle).^2));
    end

    % Optional smoothing
    if exist('sigma', 'var') && ~isempty(sigma)
        % Gaussian smoothing
        half_width   = ceil(3 * sigma);
        x            = -half_width:half_width;
        gauss_kernel = exp(-x.^2 / (2 * sigma^2));
        gauss_kernel = gauss_kernel / sum(gauss_kernel);

        % Circular convolution
        S_theta = conv([S_theta(end - half_width + 1:end), S_theta, S_theta(1:half_width)], ...
                       gauss_kernel, 'same');
        S_theta = S_theta(half_width + 1:end - half_width);
    elseif exist('windowSize', 'var') && ~isempty(windowSize)
        % Moving average smoothing
        pad    = floor(windowSize / 2);
        kernel = ones(1, windowSize) / windowSize;
        S_theta = conv([S_theta(end - pad + 1:end), S_theta, S_theta(1:pad)], kernel, 'same');
        S_theta = S_theta(pad + 1:end - pad);
    end
end