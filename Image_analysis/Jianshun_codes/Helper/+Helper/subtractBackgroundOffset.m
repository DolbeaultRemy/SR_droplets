function ret = subtractBackgroundOffset(img, fraction)
%% subtractBackgroundOffset
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Remove the background from the image.
%
% Notes:
%   Optional notes, references.
    
    x_fraction = fraction(1);
    y_fraction = fraction(2);
    offset = Helper.getBkgOffsetFromCorners(img, x_fraction, y_fraction);
    ret    = img - offset;
end