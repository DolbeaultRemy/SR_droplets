function [optrefimages] = removeFringesInImage(absimages, refimages, bgmask)
%% removefringesInImage
% Authors:      Shannon Whitlock, Caspar Ockeloen
% Date:         May 2009 (Last revision: 11 August 2010)
% Version:      Unknown
%
% Description:
%   Fringe removal and noise reduction from absorption images.
%   Creates an optimal reference image for each absorption image in a set as 
%   a linear combination of reference images, with coefficients chosen to 
%   minimize the least-squares residuals between each absorption image and 
%   the optimal reference image. The coefficients are obtained by solving a 
%   linear set of equations using matrix inverse by LU decomposition.
%
%   Application of the algorithm is described in C. F. Ockeloen et al, Improved 
%   detection of small atom numbers through image processing, arXiv:1007.2136 (2010).
%
% Inputs:
%    REQUIRED: 
%    absimages    - Absorption image data, 
%                   typically 16 bit grayscale images
%    refimages    - Raw reference image data
%       absimages and refimages are both cell arrays containing
%       2D array data. The number of refimages can differ from the 
%       number of absimages.
%
%    OPTIONAL:
%    bgmask       - Array specifying background region used, 
%                   1=background, 0=data. Defaults to all ones.
% Outputs:
%    optrefimages - Cell array of optimal reference images, 
%                   equal in size to absimages.
%
% Notes:
%   Syntax        - [optrefimages] = removefringesInImage(absimages,refimages,bgmask);
%   Dependencies  - none
%   Reference     - C. F. Ockeloen, A. F. Tauschinsky, R. J. C. Spreeuw, and
%                   S. Whitlock, Improved detection of small atom numbers through 
%                   image processing, arXiv:1007.2136
    
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