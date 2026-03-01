function results = runFitProgressViewer(od_imgs, model, quantities, plotConfig, extraParams)
%% runFitProgressViewer
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Generic batch fit viewer for 2D OD images using any model that
%   implements guess() and fit(). Dynamically updates image, fit, residual,
%   and optional bottom-row metrics (atom number, condensate fraction, temperature)
%
% Inputs:
%   od_imgs      : cell array of 2D images to fit
%   model        : object implementing guess() and fit()
%   quantities   : cell array of field names to compute/display
%   plotConfig   : struct controlling figure appearance and bottom row plots
%       .fontName, .fontSize, .fontWeight, .colormapName, .scatterLineSpec
%       .bottomRowLabels   : corresponding y-axis labels
%       .bottomRowUnits    : scaling for display
%       .bottomRowTitles   : optional subplot titles (default = labels)
%   extraParams  : struct of optional parameters for the model (e.g., ToF)
%
% Notes:
%   Optional notes, references.

    arguments
        od_imgs cell
        model
        quantities cell
        plotConfig struct = struct()
        extraParams struct = struct()
    end

    %% --- Robust defaults + user override ---

    config = struct();  % start empty
    
    % --- List of default values in a cell array ---
    defaultValues = { ...
        'fontName',        'Arial'; ...
        'fontSize',        16; ...
        'fontWeight',      'bold'; ...
        'colormapName',    'turbo'; ...
        'scatterLineSpec', '-o'; ...
        'bottomRowFields', quantities; ...
        'bottomRowLabels', {'Atom Number','# Condensed','Temp'}; ...
        'bottomRowUnits',  [1,1,1] ...
    };
    
    % --- Assign defaults field by field ---
    for k = 1:size(defaultValues,1)
        field = defaultValues{k,1};
        value = defaultValues{k,2};
        config.(field) = value;
    end
    
    % --- Merge user config if provided ---
    if exist('plotConfig','var') && isstruct(plotConfig) && ~isempty(plotConfig)
        % If plotConfig is a struct array, take first element
        if numel(plotConfig) > 1
            warning('plotConfig is a struct array. Using the first element.');
            plotConfig = plotConfig(1);
        end
        plotFields = fieldnames(plotConfig);
        for k = 1:numel(plotFields)
            f = plotFields{k};
            v = plotConfig.(f);
            % Only override if user provided a non-empty value
            if ~isempty(v)
                config.(f) = v;
            end
        end
    end
    
    % --- Safety fallback ---
    if ~isfield(config,'bottomRowFields') || isempty(config.bottomRowFields)
        config.bottomRowFields = quantities;
    end
    
    numImages = numel(od_imgs);
    fprintf('\n[INFO] Starting processing of %d images...\n', numImages);

    %% --- Preallocate results struct ---
    results = repmat(struct('imageIndex',[],'fitResult',[],'gof',[], ...
        'params',[],'fitData',[],'residuals',[],'rsquare',[],'status','Not processed'), ...
        numImages,1);

    for i = 1:numImages
        results(i).imageIndex = i;
    end

    %% --- Create or reuse figure ---
    hFig = findobj('Type','figure','Tag','FitProgressViewer');
    if isempty(hFig)
        hFig = figure('Position',[100,100,1450,850], ...
                      'NumberTitle','off', ...
                      'Name','Fit Progress Viewer', ...
                      'Tag','FitProgressViewer');
    else
        figure(hFig); clf;
    end
    t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    %% --- Pre-create image/fit/residual axes ---
    [axOriginal,hOriginal] = createImageAxis(nexttile(1), 'OD Image');
    [axFit,hFit]       = createImageAxis(nexttile(2), 'Fit');
    [axResidual,hResidual] = createImageAxis(nexttile(3), 'Residual');

    %% --- Pre-create bottom-row axes ---
    axBottom = gobjects(numel(config.bottomRowFields),1);
    scatterBottom = gobjects(numel(config.bottomRowFields),1);
    for k = 1:numel(config.bottomRowFields)
        axBottom(k) = nexttile(3+k);
        scatterBottom(k) = plot(axBottom(k), nan, nan, config.scatterLineSpec);
        hold(axBottom(k),'on'); grid(axBottom(k),'on');
        xlabel(axBottom(k),'Image Index','FontName',config.fontName);
        ylabel(axBottom(k),config.bottomRowLabels{k},'FontName',config.fontName);

        % Use bottomRowTitles if provided, else fallback to labels
        if isfield(config,'bottomRowTitles') && numel(config.bottomRowTitles) >= k
            titleStr = config.bottomRowTitles{k};
        else
            titleStr = config.bottomRowLabels{k};
        end
        title(axBottom(k), titleStr, 'FontName',config.fontName, ...
            'FontSize',config.fontSize, 'FontWeight',config.fontWeight);
    end

    %% --- Apply consistent font formatting ---
    allAxes = [axOriginal, axFit, axResidual, axBottom(:)'];
    for ax = allAxes, set(ax,'FontSize',config.fontSize,'FontName',config.fontName); end

    %% --- Main batch loop ---
    for i = 1:numImages
        currentImg = od_imgs{i};
        if isempty(currentImg) || ~isnumeric(currentImg) || all(isnan(currentImg(:)))
            warning('Image %d empty or invalid. Skipping',i);
            results(i).status = 'Invalid image'; continue;
        end

        [ny,nx] = size(currentImg); x = 1:nx; y = 1:ny;

        %% --- Model guess and fit ---
        params = model.guess(currentImg,x,y);

        if isempty(fieldnames(extraParams))
            [fitResult,gof] = model.fit(currentImg,x,y,'params',params);
        else
            args = reshape([fieldnames(extraParams)'; struct2cell(extraParams)'],1,[]);
            [fitResult,gof] = model.fit(currentImg,x,y,'params',params,args{:});
        end

        [X,Y] = meshgrid(x,y);
        xyData = [X(:),Y(:)];
        fitData = reshape(fitResult(xyData), size(currentImg));
        residuals = currentImg - fitData;

        %% --- Store results ---
        results(i).fitResult = fitResult;
        results(i).gof       = gof;
        results(i).params    = params;
        results(i).fitData   = fitData;
        results(i).residuals = residuals;
        results(i).rsquare   = gof.rsquare;
        results(i).status    = 'Success';

        %% --- Compute requested bottom-row fields only ---
        for k = 1:numel(config.bottomRowFields)
            fieldName = config.bottomRowFields{k};
            switch fieldName
                case 'atom_number'
                    if ismethod(model,'return_atom_number')
                        atomStruct = model.return_atom_number(X,Y,false);
                        results(i).atom_number = atomStruct.N_bec;
                    end
                case 'condensate_fraction'
                    if isprop(model,'cond_frac')
                        results(i).condensate_fraction = model.cond_frac;
                    end
                case 'temperature'
                    if ismethod(model,'return_temperature') && isfield(extraParams,'ToF')
                        results(i).temperature = model.return_temperature(extraParams.ToF,[],false);
                    end
            end
        end

        %% --- Update plots dynamically ---
        updatePlots(currentImg, fitData, residuals, i);
    end

    %% --- Display mean ± SEM after full batch loop (raw values) ---
    for k = 1:numel(config.bottomRowFields)
        fieldName = config.bottomRowFields{k};
        
        % Collect raw values (no scaling)
        vals = nan(numImages,1);
        for i = 1:numImages
            if isfield(results(i), fieldName)
                vals(i) = results(i).(fieldName);  % RAW, unscaled
            end
        end
        
        % Only keep valid entries
        validVals = vals(~isnan(vals));
        nVals = numel(validVals);
    
        meanVal = mean(validVals,'omitnan');
        
        if nVals >= 20
            % Use SEM
            semVal = std(validVals,'omitnan')/sqrt(nVals);
            str = sprintf('Mean ± SEM: %.2e ± %.2e', meanVal, semVal);
        else
            % Use SD
            sdVal = std(validVals,'omitnan');
            str = sprintf('Mean ± SD: %.2e ± %.2e', meanVal, sdVal);
        end
        
        % Place in bottom-right corner using normalized axes coordinates
        ax = axBottom(k);
        text(ax, 0.98, 0.02, str, ...
             'Units','normalized', ...
             'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
             'FontName', config.fontName, ...
             'FontSize', config.fontSize-2, ...
             'FontWeight', config.fontWeight, ...
             'BackgroundColor', 'w', ...          % white box
             'Margin', 4, ...                     % padding inside box
             'EdgeColor', 'k');                   % black border
    end
    
    fprintf('\n[INFO] Processing complete.\n');

    %% --- Nested functions ---
    function [ax,hImg] = createImageAxis(parentTile, titleStr)
        ax = parentTile;
        hImg = imagesc(ax, nan); axis(ax,'equal','tight');
        colormap(ax, config.colormapName); colorbar(ax);
        title(ax, titleStr, 'FontName', config.fontName, ...
            'FontSize', config.fontSize, 'FontWeight', config.fontWeight);
    end

    function updatePlots(img, fitImg, residImg, idx)
        % Update main image axes
        set(hOriginal, 'CData', img);
        set(hFit,      'CData', fitImg);
        set(hResidual, 'CData', residImg);

        % Update bottom-row plots
        for k = 1:numel(config.bottomRowFields)
            fieldName = config.bottomRowFields{k};
            val = NaN;
            if isfield(results(idx), fieldName)
                val = results(idx).(fieldName) * config.bottomRowUnits(k);
            end
            if idx == 1
                scatterBottom(k).XData = idx;
                scatterBottom(k).YData = val;
            else
                scatterBottom(k).XData = [scatterBottom(k).XData, idx];
                scatterBottom(k).YData = [scatterBottom(k).YData, val];
            end
        end
        drawnow;
    end
end
