function [theta_vals, g2_angular] = computeAngularCorrelation(g2_matrix, dx_phys, dy_phys, r_min, r_max, num_bins)
%% computeAngularCorrelation
% Author:       Karthik
% Date:         2025-09-14
% Version:      1.0
%
% Description:
%   Extracts angular profile of g2(Δx,Δy) along a radial band [r_min, r_max]
%   from 0 to 180 degrees.
%
% Inputs:
%   g2_matrix : 2D autocorrelation matrix
%   dx_phys   : x-axis shifts in µm
%   dy_phys   : y-axis shifts in µm
%   r_min     : minimum radial distance (µm)
%   r_max     : maximum radial distance (µm)
%   num_bins  : number of angular bins
%
% Outputs:
%   theta_vals  : angular positions [radians]
%   g2_angular  : angular profile of g2
%
% Notes:
%   Optional notes, references.

    [X, Y] = meshgrid(dx_phys, dy_phys);
    R = sqrt(X.^2 + Y.^2);
    Theta = atan2(Y, X);  % [-pi, pi]

    % Restrict to the radial band
    radial_mask = (R >= r_min) & (R <= r_max);

    % Angular bins from 0 to pi (0°–180°)
    theta_vals = linspace(0, pi, num_bins);
    g2_angular = zeros(1, num_bins);

    for i = 1:num_bins
        angle_start = theta_vals(i) - (theta_vals(2)-theta_vals(1))/2;
        angle_end   = theta_vals(i) + (theta_vals(2)-theta_vals(1))/2;

        % Handle wrap-around at 0/pi
        if angle_start < 0
            angle_mask = (Theta >= 0 & Theta < angle_end) | (Theta >= (pi+angle_start) & Theta <= pi);
        elseif angle_end > pi
            angle_mask = (Theta >= angle_start & Theta <= pi) | (Theta >= 0 & Theta < (angle_end - pi));
        else
            angle_mask = (Theta >= angle_start) & (Theta < angle_end);
        end

        % Combine with radial mask
        bin_mask = radial_mask & angle_mask;

        % Sum or average within this angular bin
        g2_angular(i) = mean(g2_matrix(bin_mask), 'omitnan');
    end
end
