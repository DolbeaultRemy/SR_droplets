function results = conductPCA(od_imgs)
%% computePCAfromImages
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Performs PCA on optical density images and returns results in a struct.
%
% Inputs:
%   od_imgs - cell array of OD images
%
% Outputs:
%   pcaResults - struct containing PCA outputs:
%       .coeff    - PCA coefficients (principal components)
%       .score    - PCA scores for each image
%       .explained - variance explained by each PC
%       .Nx, .Ny   - dimensions of individual images
%
% Notes:
%   Optional notes, references.

    allImgs3D = cat(3, od_imgs{:});
    [Nx, Ny]  = size(allImgs3D(:,:,1));
    Xall      = reshape(allImgs3D, [], numel(od_imgs))';
    [coeff, score, ~, ~, explained] = pca(Xall);

    results = struct( ...
        'coeff', coeff, ...
        'score', score, ...
        'explained', explained, ...
        'Nx', Nx, ...
        'Ny', Ny ...
    );
end
