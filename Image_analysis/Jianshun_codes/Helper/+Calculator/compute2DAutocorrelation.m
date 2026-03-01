function [g2, dx_phys, dy_phys] = compute2DAutocorrelation(img, max_shift, pixel_size, magnification)
%% compute2DAutocorrelation
% Author:       Karthik
% Date:         2025-09-14
% Version:      1.0
%
% Description:
%   Computes 2D autocorrelation g2(Δx,Δy) with shifts in physical units [µm]
%
% Inputs:
%   img          : 2D image (OD)
%   max_shift_um : maximum shift in micrometers
%   pixel_size   : camera pixel size in meters
%   magnification: imaging magnification
%
% Outputs:
%   g2           : normalized 2D autocorrelation
%   dx_phys      : x-axis shifts in µm
%   dy_phys      : y-axis shifts in µm
%
% Notes:
%   Optional notes, references.

    [M, N]       = size(img);
    img_mean_sub = img - mean(img(:));
    I2_mean      = mean(img_mean_sub(:).^2);
    
    % Convert max_shift from µm to pixels
    dx_max_px = round(max_shift / (pixel_size/magnification * 1e6));
    
    g2 = zeros(2*dx_max_px+1, 2*dx_max_px+1);
    
    for dx = -dx_max_px:dx_max_px
        for dy = -dx_max_px:dx_max_px
            % overlapping region in pixels
            x1  = max(1,1+dx); x2 = min(M,M+dx);
            y1  = max(1,1+dy); y2 = min(N,N+dy);
            x1s = max(1,1-dx); x2s = min(M,M-dx);
            y1s = max(1,1-dy); y2s = min(N,N-dy);
            
            overlap = img_mean_sub(x1:x2, y1:y2) .* img_mean_sub(x1s:x2s, y1s:y2s);
            g2(dx+dx_max_px+1, dy+dx_max_px+1) = mean(overlap(:)) / I2_mean;
        end
    end
    
    % Return physical shift axes in µm
    dx_phys = (-dx_max_px:dx_max_px) * (pixel_size/magnification * 1e6);
    dy_phys = dx_phys;
end
