function ret = cropODImage(img, center, span)
%% cropODImage
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Crop the image according to the region of interest (ROI).
%
% Notes:
%   Optional notes, references.

    x_start = floor(center(1) - span(1) / 2);
    x_end   = floor(center(1) + span(1) / 2);
    y_start = floor(center(2) - span(2) / 2);
    y_end   = floor(center(2) + span(2) / 2);

    ret = img(y_start:y_end, x_start:x_end);
end