function saveFigure(fig, varargin)
%% saveFigure
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Saves a MATLAB figure as a .fig file in a specified directory.
%
% Inputs:
%   fig             - Figure handle to save
%
% Optional Parameters:
%   'SaveFileName'  - Name of the file (default: 'figure.fig')
%   'SaveDirectory' - Directory to save into (default: current working directory)
%   'SkipSaveFigures' - If true, skips saving (default: false)
%
% Notes:
% Usage -
%   saveFigure(fig)
%   saveFigure(fig, 'SaveFileName', 'myplot.fig', 'SaveDirectory', 'results', 'SkipSaveFigures', false)
%
% Example -
%   fig = figure;
%   plot(1:10, rand(1,10));
%   saveFigure(fig, 'SaveFileName', 'test.fig', 'SaveDirectory', 'plots');

    % --- Defaults ---
    p = inputParser;
    addParameter(p, 'SaveFileName', 'figure.fig');
    addParameter(p, 'SaveDirectory', pwd);
    addParameter(p, 'SkipSaveFigures', false);
    parse(p, varargin{:});
    opts = p.Results;

    if opts.SkipSaveFigures
        return; % Do nothing
    end

    % --- Ensure directory exists ---
    if ~exist(opts.SaveDirectory, 'dir')
        mkdir(opts.SaveDirectory);
    end

    % --- Ensure .fig extension ---
    [~, name, ext] = fileparts(opts.SaveFileName);
    if isempty(ext)
        ext = '.fig';
    elseif ~strcmpi(ext, '.fig')
        warning('Overriding extension to .fig (was %s).', ext);
        ext = '.fig';
    end

    saveFullPath = fullfile(opts.SaveDirectory, [name ext]);
    savefig(fig, saveFullPath);
    fprintf('Figure saved as MATLAB .fig: %s\n', saveFullPath);
end
