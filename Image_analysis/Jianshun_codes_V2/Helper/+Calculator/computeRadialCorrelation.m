function [r_vals, g2_profile] = computeRadialCorrelation(g2_matrix, dx_phys, dy_phys, theta0)
%% computeRadialCorrelation
% Author:       Karthik
% Date:         2025-09-14
% Version:      1.0
%
% Description:
%   Extracts the g2 profile along a specific radial direction (line) at angle theta0
%
% Inputs:
%   g2_matrix   : 2D autocorrelation matrix
%   dx_phys     : x-axis shifts in µm
%   dy_phys     : y-axis shifts in µm
%   theta0_deg  : angle of the radial direction (degrees)
%
% Outputs:
%   r_vals      : radial distances along the line (µm)
%   g2_profile  : g2 values along that line
%
% Notes:
%   Optional notes, references.

    % Create meshgrid of physical shifts
    [X, Y] = meshgrid(dx_phys, dy_phys);

    % Compute radial distance along line
    R = sqrt(X.^2 + Y.^2);
    
    % Compute angle at each point
    Theta = atan2(Y, X);

    % Mask for points along the chosen line
    % Use a small tolerance to pick points close to the line
    tol = 0.5 * mean(diff(dx_phys));  % half pixel tolerance
    line_mask = abs(Theta - theta0) < tol;

    % Extract g2 values along the line
    r_vals = R(line_mask);
    g2_profile = g2_matrix(line_mask);

    % Sort by radial distance
    [r_vals, sort_idx] = sort(r_vals);
    g2_profile = g2_profile(sort_idx);
end
