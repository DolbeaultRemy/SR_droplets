function [od_imgs,Date] = run_data_remy(folderPath, run, cam, groupList, span, center, fraction ,shotwindow, removeFringes)
    

    %% input example:
    % folderPath   = "\\DyLabNAS\Data\TwoDGas\2025\12\02\\";
    % run          = "0000";
    % groupList    = ["/images/MOT_3D_Camera/in_situ_absorption", "/images/ODT_1_Axis_Camera/in_situ_absorption", "/images/ODT_2_Axis_Camera/in_situ_absorption", "/images/Horizontal_Axis_Camera/in_situ_absorption", "/images/Vertical_Axis_Camera/in_situ_absorption"];
    % cam          = 5; %% (vertical)
    % span         = [150, 150]; %% ROI stuff (I larger but if zou are mainly doing in-situ... )
    % center       = [1420,2050]; %% ROI stuff (For recent times should be around here)
    % fraction     = [0.1, 0.1]; %% Chooses fraction in corner
    % shotwindow   = 1 %% amountof shots in each run
    % removeFringes = false;

    folderPath     = strcat(folderPath, run);    
    filePattern = fullfile(folderPath, '*.h5');
    files       = dir(filePattern);
    k = 1;
    textprogressbar('Reading Data: ');
    Date        = files(1).date;
    for j = 1:shotwindow
        baseFileName = files(j).name;
        fullFileName = fullfile(files(j).folder, baseFileName);

        % fprintf(1, 'Now reading %s\n', fullFileName);
        
        atm_img  = double(imrotate(h5read(fullFileName, append(groupList(cam), "/atoms")), 0)); % im2double rescales values to between [0, 1], use double instead
        bkg_img  = double(imrotate(h5read(fullFileName, append(groupList(cam), "/background")), 0));
        dark_img = double(imrotate(h5read(fullFileName, append(groupList(cam), "/dark")), 0));

        refimages(:,:,k)  = subtractBackgroundOffset(cropODImage(bkg_img, center, span), fraction)';
        absimages(:,:,k)  = subtractBackgroundOffset(cropODImage(calculateODImage(atm_img, bkg_img, dark_img), center, span), fraction)';
        
        textprogressbar(k/shotwindow*100)
        k = k + 1;
    end
    textprogressbar(' Done');
    if removeFringes
        fprintf(1, '--Removing Fringes--\n');
        optrefimages               = removefringesInImage(absimages, refimages);
        absimages_fringe_removed   = absimages(:, :, :) - optrefimages(:, :, :);

        nimgs                      = size(absimages_fringe_removed,3);
        od_imgs                    = cell(1, nimgs);
        for i = 1:nimgs
            od_imgs{i}             = absimages_fringe_removed(:, :, i);
        end
    else
        nimgs                      = shotwindow;
        od_imgs                    = cell(1, nimgs);
        for i = 1:nimgs
            od_imgs{i}             = smoothdata(absimages(:, :, i));
        end
    end
end


%%% Helper functions

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
