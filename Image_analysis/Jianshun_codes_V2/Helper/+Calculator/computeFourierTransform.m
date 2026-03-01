function [IMGFFT, IMGPR] = computeFourierTransform(I, skipPreprocessing, skipMasking, skipIntensityThresholding, skipBinarization)
%% computeFourierTransform
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Computes the 2D Fourier power spectrum of binarized and enhanced lattice image features, with optional central mask.
% Inputs:
%   I           - Grayscale or RGB image matrix
%
% Output:
%   F_mag       - 2D Fourier power spectrum (shifted)
%
% Notes:
%   Optional notes, references.
    
    if ~skipPreprocessing
        % Preprocessing: Denoise
        filtered             = imgaussfilt(I, 10);
        IMGPR                = I - filtered; % adjust sigma as needed
    else
        IMGPR                = I;
    end
    
    if ~skipMasking
        [rows, cols]     = size(IMGPR);
        [X, Y]           = meshgrid(1:cols, 1:rows);
        % Elliptical mask parameters
        cx               = cols / 2;
        cy               = rows / 2;
        
        % Shifted coordinates
        x                = X - cx;
        y                = Y - cy;
        
        % Ellipse semi-axes
        rx               = 0.4 * cols;
        ry               = 0.2 * rows;
        
        % Rotation angle in degrees -> radians
        theta_deg        = 30;        % Adjust as needed
        theta            = deg2rad(theta_deg);
        
        % Rotated ellipse equation
        cos_t            = cos(theta);
        sin_t            = sin(theta);
        
        x_rot            = (x * cos_t + y * sin_t);
        y_rot            = (-x * sin_t + y * cos_t);
        
        ellipseMask      = (x_rot.^2) / rx^2 + (y_rot.^2) / ry^2 <= 1;
        
        % Apply cutout mask
        IMGPR            = IMGPR .* ellipseMask;
    end

    if ~skipIntensityThresholding
        % Apply global intensity threshold mask
        intensity_thresh = 0.20; 
        intensity_mask   = IMGPR > intensity_thresh;
        IMGPR            = IMGPR .* intensity_mask;
    end

    if ~skipBinarization
        % Adaptive binarization and cleanup
        IMGPR            = imbinarize(IMGPR, 'adaptive', 'Sensitivity', 0.0);
        IMGPR            = imdilate(IMGPR, strel('disk', 2));
        IMGPR            = imerode(IMGPR, strel('disk', 1));
        IMGPR            = imfill(IMGPR, 'holes');
        F                = fft2(double(IMGPR)); % Compute 2D Fourier Transform
        IMGFFT           = abs(fftshift(F))';  % Shift zero frequency to center
    else
        F                = fft2(double(IMGPR)); % Compute 2D Fourier Transform
        IMGFFT           = abs(fftshift(F))';  % Shift zero frequency to center
    end
end