% folderPath = "C:\Users\Joschka\Desktop\Remy\Codes\SR_droplets\Image_analysis\Jianshun_codes\Data_analysed\" + year + "\" + month + "\" + day + "\" + shot + "\";
folderPath = "Y:\StructuralPhaseTransition\2026\02\11\";
ROIcenter = [1500 1950];
measurementName = "Stripes_90°";
run = "0001";
cam = 5;
groupList    = ["/images/MOT_3D_Camera/in_situ_absorption", "/images/ODT_1_Axis_Camera/in_situ_absorption", "/images/ODT_2_Axis_Camera/in_situ_absorption", "/images/Horizontal_Axis_Camera/in_situ_absorption", "/images/Vertical_Axis_Camera/in_situ_absorption"];
span         = [400, 400];
center       = [1500, 1950];
fraction     = [0.1, 0.1];
shotwindow   = 5;%380;
removeFringes = false;

[od_imgs,Date] = run_data_remy(folderPath, run, cam, groupList, span, center, fraction ,shotwindow, removeFringes);

detectionParams = setDetectionParameters();

folderPath     = strcat(folderPath, run);    
filePattern = fullfile(folderPath, '*.h5');
files       = dir(filePattern);


%%% Detect structures and create masked images %%%
% Initialize cell arrays to store results
masks_detected_structures = {};
xStart_mask = {};
yStart_mask = {};
od_imgs_crop = {};
masked_images_atoms = {};

% Process all images
for idx = 1:length(od_imgs)
    % Detect structures (your existing code)
    [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart, binaryMask, CC] = ...
        detectStructure(od_imgs{1, idx}, detectionParams, folderPath, files(idx).name);
    
    % Store detection results
    masks_detected_structures{end+1} = binaryMask;
    xStart_mask{end+1} = xStart;
    yStart_mask{end+1} = yStart;
    % 
    % % Create full-sized binary mask matching original image dimensions
    % [orig_height, orig_width] = size(od_imgs{idx});
    % full_binary_mask = false(orig_height, orig_width);
    % 
    % % Place the cropped binary mask at the correct position
    % y_end = min(yStart + size(binaryMask, 1) - 1, orig_height);
    % x_end = min(xStart + size(binaryMask, 2) - 1, orig_width);
    % full_binary_mask(yStart:y_end, xStart:x_end) = binaryMask(1:(y_end-yStart+1), 1:(x_end-xStart+1));
    % 
end

% % Save masked_images to HDF5 file with custom name %%%
% 
% % Set your desired filename here
% customFilename = 'my_masked_images.h5';  % CHANGE THIS TO YOUR DESIRED NAME
% 
% % Create full path
% outputFilename = fullfile(customFilename);
% 
% % Convert to character vector if needed
% outputFilename = char(outputFilename);
% 
% % Delete existing file if it exists
% if exist(outputFilename, 'file') == 2
%     delete(outputFilename);
%     fprintf('Deleted existing file: %s\n', outputFilename);
% end
% 
% % Save each masked image as a separate dataset
% for idx = 1:length(masked_images)
%     % Create dataset name
%     datasetName = sprintf('/image_%04d', idx);
%     datasetName = char(datasetName);
% 
%     % Save the image
%     h5create(outputFilename, datasetName, size(masked_images{idx}), 'Datatype', 'double');
%     h5write(outputFilename, datasetName, masked_images{idx});
% 
%     % Optional: save minimal metadata
%     h5writeatt(outputFilename, datasetName, 'original_index', idx);
% end
% 
% % Save basic info about the dataset
% h5writeatt(outputFilename, '/', 'num_images', length(masked_images));
% h5writeatt(outputFilename, '/', 'description', 'Masked images with non-detected pixels set to 0');
% 
% fprintf('Saved %d masked images to: %s\n', length(masked_images), outputFilename);







function params = setDetectionParameters()
% Set parameters for structure detection.
    params.backgroundDiskFraction    = 0.1250;
    params.boundingBoxPadding        = 35;
    params.dogGaussianSmallSigma     = 0.5;
    params.dogGaussianLargeSigma     = 4;
    params.adaptiveSensitivity       = 0.3;
    params.adaptiveNeighborhoodSize  = 13;
    params.minPeakFraction           = 0.2;
    params.minimumPatchArea          = 20;
    params.shapeMinArea              = 20;
    params.shapeCloseRadius          = 3;
    params.shapeFillHoles            = false;
    params.intensityThreshFraction   = 0.4499;
    params.edgeSigma                 = 1.1749;
    params.edgeThresholdLow          = 0.3383;
    params.edgeThresholdHigh         = 0.6412;
    params.pixelSize                 = 5.86e-6;
    params.magnification             = 23.94;
end

%%% Helper functions



function ret = subtractBackgroundOffset(img, fraction)
% Remove the background from the image.
% :param dataArray: The image
% :type dataArray: xarray DataArray
% :param x_fraction: The fraction of the pixels used in x axis
% :type x_fraction: float
% :param y_fraction: The fraction of the pixels used in y axis
% :type y_fraction: float
% :return: The image after removing background
% :rtype: xarray DataArray

x_fraction = fraction(1);
y_fraction = fraction(2);
offset = getBkgOffsetFromCorners(img, x_fraction, y_fraction);
ret    = img - offset;
end

function ret = getBkgOffsetFromCorners(img, x_fraction, y_fraction)
% image must be a 2D numerical array
[dim1, dim2] = size(img);

s1 = img(1:round(dim1 * y_fraction), 1:round(dim2 * x_fraction));
s2 = img(1:round(dim1 * y_fraction), round(dim2 - dim2 * x_fraction):dim2);
s3 = img(round(dim1 - dim1 * y_fraction):dim1, 1:round(dim2 * x_fraction));
s4 = img(round(dim1 - dim1 * y_fraction):dim1, round(dim2 - dim2 * x_fraction):dim2);

ret = mean([mean(s1(:)), mean(s2(:)), mean(s3(:)), mean(s4(:))]);
end

function ret = calculateODImage(imageAtom, imageBackground, imageDark)
% Calculate the OD image for absorption imaging.
% :param imageAtom: The image with atoms
% :type imageAtom: numpy array
% :param imageBackground: The image without atoms
% :type imageBackground: numpy array
% :param imageDark: The image without light
% :type imageDark: numpy array
% :return: The OD images
% :rtype: numpy array

numerator   = imageBackground - imageDark;
denominator = imageAtom - imageDark;

numerator(numerator == 0)     = 1;
denominator(denominator == 0) = 1;

% ret = -log(double(abs(denominator ./ numerator)));
exposureTime = 5e-6;
imageOD = double(abs(denominator ./ numerator));
ret = -log(imageOD) + (numerator - denominator) ./ (7000 * (exposureTime / 5e-6));
if numel(ret) == 1
    ret = ret(1);
end
end

function [optrefimages] = removefringesInImage(absimages, refimages, bgmask)
% removefringesInImage - Fringe removal and noise reduction from absorption images.
% Creates an optimal reference image for each absorption image in a set as
% a linear combination of reference images, with coefficients chosen to
% minimize the least-squares residuals between each absorption image and
% the optimal reference image. The coefficients are obtained by solving a
% linear set of equations using matrix inverse by LU decomposition.
%
% Application of the algorithm is described in C. F. Ockeloen et al, Improved
% detection of small atom numbers through image processing, arXiv:1007.2136 (2010).
%
% Syntax:
%    [optrefimages] = removefringesInImage(absimages,refimages,bgmask);
%
% Required inputs:
%    absimages    - Absorption image data,
%                   typically 16 bit grayscale images
%    refimages    - Raw reference image data
%       absimages and refimages are both cell arrays containing
%       2D array data. The number of refimages can differ from the
%       number of absimages.
%
% Optional inputs:
%    bgmask       - Array specifying background region used,
%                   1=background, 0=data. Defaults to all ones.
% Outputs:
%    optrefimages - Cell array of optimal reference images,
%                   equal in size to absimages.
%

% Dependencies: none
%
% Authors: Shannon Whitlock, Caspar Ockeloen
% Reference: C. F. Ockeloen, A. F. Tauschinsky, R. J. C. Spreeuw, and
%            S. Whitlock, Improved detection of small atom numbers through
%            image processing, arXiv:1007.2136
% Email:
% May 2009; Last revision: 11 August 2010

% Process inputs

% Set variables, and flatten absorption and reference images
nimgs  = size(absimages,3);
nimgsR = size(refimages,3);
xdim = size(absimages(:,:,1),2);
ydim = size(absimages(:,:,1),1);

R = single(reshape(refimages,xdim*ydim,nimgsR));
A = single(reshape(absimages,xdim*ydim,nimgs));
optrefimages=zeros(size(absimages)); % preallocate

if not(exist('bgmask','var')); bgmask=ones(ydim,xdim); end
k = find(bgmask(:)==1);  % Index k specifying background region

% Ensure there are no duplicate reference images
% R=unique(R','rows')'; % comment this line if you run out of memory

% Decompose B = R*R' using singular value or LU decomposition
[L,U,p] = lu(R(k,:)'*R(k,:),'vector');       % LU decomposition

for j=1:nimgs
    b=R(k,:)'*A(k,j);
    % Obtain coefficients c which minimise least-square residuals
    lower.LT = true; upper.UT = true;
    c = linsolve(U,linsolve(L,b(p,:),lower),upper);

    % Compute optimised reference image
    optrefimages(:,:,j)=reshape(R*c,[ydim xdim]);
end
end