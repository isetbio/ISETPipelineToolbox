function buildMosaicMontage(pr, cnv, stage, varargin)
% Synopsis:
%    Build render matrix if desired/needed
%
% Description:
%    [UPDATE]. Run this function if you would like to rebuild a new mosaic and
%    render matrix.  This also gets run if there is no cached file corresponding
%    to the desired parameters. Once built, this file can be loaded from cache
%    for quicker running.
%
% See also: aoStimRecon, aoStimReconRunMany, buildBaseMosaic,
%           buildRenderStruct.

% History:
%   08/16/22  chr  Made it a callable function
%   08/25/22  chr  Included portion for dichromacy
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options
%   03/13/24  chr  Split from buildRenderStruct

%% Set the stage and build base mosaic
%
% Unpack the variables that change dependent on whether building a forward
% mosaic or recon mosaic. useDisplay in reference to whether a display
% needs to actually be loaded to achieve desired function (true for render
% struct, false for montage)
useDisplay = false;
st = unpackStage(pr, cnv, stage);
[theConeMosaic, ~] = buildBaseMosaic(pr, cnv, st, useDisplay);

%% Axes Naming
%
% Adjust how variables are presented to make montage prettier
for w = 1:length(pr.focalRegionList)
    switch pr.focalRegionList(w)
        case "center"
            regionAxesNames(w) = "Center";
        case "nearSurround"
            regionAxesNames(w) = "Near";
        case "distantSurround"
            regionAxesNames(w) = "Distant";
        case "multiple"
            regionAxesNames(w) = "Multiple";
        case "global"
            regionAxesNames(w) = "Global";
    end
end

switch stage
    case "forward"
        stageTitle = "Forward";
    case "recon"
        stageTitle = "Recon";
end

%% Build mosaic montage based on edited mosaic
%
% Establish dimensions of the full montage
allRows = length(pr.stimSizeDegsList) * length(pr.focalRegionList);
allColms = length(pr.focalPropLList);
viewBounds = false; %%%%%%%%%%%%%

allConeMosaics = cell(allRows, allColms);
allMosaicsConeInfo = cell(allRows, allColms);
allRenderNames = cell(allRows, allColms);

if pr.useCustomMosaic
    for h = 1:length(pr.focalVariantList)

        theFig = figure;
        t = tiledlayout(allRows, allColms, 'TileSpacing','none');


        for i = 1:length(pr.stimSizeDegsList)
            for j = 1:length(pr.focalRegionList)
                for k = 1:length(pr.focalPropLList)

                    [theConeMosaic, mosaicConeInfo] = setConeProportions(pr.focalRegionList(j), ...
                        pr.focalPropLList(k), pr.focalVariantList(h), theConeMosaic, pr.eccXDegs, pr.eccYDegs, ...
                        pr.stimSizeDegsList(i), pr.fieldSizeMinutes, pr.regionVariant, pr.propL, pr.propS);

                    rowPlace = j + (length(pr.stimSizeDegsList) * (i - 1));
                    allConeMosaics{rowPlace,k} = theConeMosaic;
                    allMosaicsConeInfo{rowPlace,k} = mosaicConeInfo;
                    allRenderNames{rowPlace,k} = fullfile(st.renderDirFull, cnv.renderName);

                    % Plot the mosaic in the montage
                    theAxes = nexttile;
                    figureHandle = theFig;
                    theConeMosaic.visualizeMosaic(figureHandle,theAxes);
                    set(gca, 'xticklabel', [], 'yticklabel', []);
                    set(gca, 'xlabel', [], 'ylabel', []);

                    hold on;

                    if viewBounds
                        % Pull region boundary info
                        xBounds = mosaicConeInfo.xBounds;
                        yBounds = mosaicConeInfo.yBounds;

                        % Superimpose the boundaries
                        for w = length(xBounds)
                            rectangle('Position', ...
                                [xBounds(w,1) yBounds(w,1) ...
                                (xBounds(w,2) - xBounds(w,1)) ...
                                (yBounds(w,2) - yBounds(w,1))], ...
                                'LineWidth', 3)
                        end
                    end

                    % Set plot figure information
                    if i == 1 && j == 1
                        title([num2str(pr.focalPropLList(k)) 'L'])
                    end

                    if k == 1
                        ylabel({regionAxesNames(j); [num2str(pr.stimSizeDegsList(i)*60) ' Stim'] }, 'FontWeight', 'bold')
                    end
                end
            end
        end

        % Spruce up the montage figure
        set(gcf, 'Position', [595 5 1361 972]);
        title(t, {[char(stageTitle) ' Mosaic Montage'] ...
            ['Variant ' char(num2str(pr.focalVariantList(h)))] ...
            ['  ']}, ...
            'FontSize', 36)
        
        % Start creating the output file name
        outputName = ['variant', int2str(pr.focalVariantList(h)), '_version'];
      
        % Track the history of any montages already run and cached, then
        % add one to the current output file. This approach avoids
        % overwriting output results. Need to update again, make it such
        % that there is a max feature and maybe a name tag associated with
        % it. 
        versionHistory = dir(fullfile(st.montageDirFull));
        currentVersion = 1;
        for vh = 1:length(versionHistory)
            storedNames = versionHistory(vh).name;
            if contains(storedNames, outputName)
                startPoint = strfind(storedNames, 'version') + length('version'); 
                endPoint = strfind(storedNames, '.') - 1;
                storedVersion = str2double(storedNames(startPoint:endPoint));
                currentVersion = storedVersion + 1;
            end
        end
        
        saveas(gcf, fullfile(st.montageDirFull, ...
            [outputName int2str(currentVersion) '.tiff']), 'tiff');

        save(fullfile(st.montageDirFull, ...
            [outputName int2str(currentVersion) '.mat']), ...
            "allConeMosaics", "allMosaicsConeInfo", "allRenderNames");

    end
else
    disp('Not building customized mosaic montage');
end

end

