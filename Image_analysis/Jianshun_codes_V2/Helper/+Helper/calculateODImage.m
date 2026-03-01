function imageOD = calculateODImage(imageAtom, imageBackground, imageDark, mode, exposureTime)
%% calculateODImage
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Calculates the optical density (OD) image for absorption imaging.
%
% Inputs:
%   imageAtom       - Image with atoms
%   imageBackground - Image without atoms
%   imageDark       - Image without light
%   mode            - 'LowIntensity' (default) or 'HighIntensity'
%   exposureTime    - Required only for 'HighIntensity' [in seconds]
%
% Output:
%   imageOD         - Computed OD image
%
% Notes:
%   Syntax - imageOD = calculateODImage(imageAtom, imageBackground, imageDark, mode, exposureTime)

    arguments
        imageAtom       (:,:) {mustBeNumeric}
        imageBackground (:,:) {mustBeNumeric}
        imageDark       (:,:) {mustBeNumeric}
        mode            char {mustBeMember(mode, {'LowIntensity', 'HighIntensity'})} = 'LowIntensity'
        exposureTime    double = NaN
    end

    % Compute numerator and denominator
    numerator   = imageBackground - imageDark;
    denominator = imageAtom - imageDark;

    % Avoid division by zero
    numerator(numerator == 0)     = 1;
    denominator(denominator == 0) = 1;

    % Calculate OD based on mode
    switch mode
        case 'LowIntensity'
            imageOD = -log(abs(denominator ./ numerator));
            
        case 'HighIntensity'
            if isnan(exposureTime)
                error('Exposure time must be provided for HighIntensity mode.');
            end
            imageOD = abs(denominator ./ numerator);
            imageOD = -log(imageOD) + (numerator - denominator) ./ (7000 * (exposureTime / 5e-6));
    end

end