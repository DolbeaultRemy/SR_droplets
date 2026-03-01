function ret = cropODImage(img, center, span)
% Crop the image according to the region of interest (ROI).
% :param dataSet: The images
% :type dataSet: xarray DataArray or DataSet
% :param center: The center of region of interest (ROI)
% :type center: tuple
% :param span: The span of region of interest (ROI)
% :type span: tuple
% :return: The cropped images
% :rtype: xarray DataArray or DataSet

x_start = floor(center(1) - span(1) / 2);
x_end   = floor(center(1) + span(1) / 2);
y_start = floor(center(2) - span(2) / 2);
y_end   = floor(center(2) + span(2) / 2);

ret = img(y_start:y_end, x_start:x_end);
end