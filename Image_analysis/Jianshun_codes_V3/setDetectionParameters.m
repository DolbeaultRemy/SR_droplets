
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
