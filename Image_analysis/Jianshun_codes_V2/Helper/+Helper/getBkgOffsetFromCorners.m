function ret = getBkgOffsetFromCorners(img, x_fraction, y_fraction)
%% getBkgOffsetFromCorners
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Brief description of the script functionality.
%
% Notes:
%   Optional notes, references.

    % image must be a 2D numerical array
    [dim1, dim2] = size(img);

    s1 = img(1:round(dim1 * y_fraction), 1:round(dim2 * x_fraction)); 
    s2 = img(1:round(dim1 * y_fraction), round(dim2 - dim2 * x_fraction):dim2);
    s3 = img(round(dim1 - dim1 * y_fraction):dim1, 1:round(dim2 * x_fraction));
    s4 = img(round(dim1 - dim1 * y_fraction):dim1, round(dim2 - dim2 * x_fraction):dim2);

    ret = mean([mean(s1(:)), mean(s2(:)), mean(s3(:)), mean(s4(:))]);
end