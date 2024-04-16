function buildSummaryPlots(pr, cnv, pp, numStim, numProp)

% If there's no such name, build this one as the first version. Otherwise
% name it as one version beyond the current existing one. Will have to
% manually delete old versions, but at least this way able to keep track of
% them for as long as desired.
%
% Also include an option to delete all previous versions in one place and
% just keep the largest one (most recent).

%% Retrieve the pertinent information cached from aoStimRecon
load(fullfile(cnv.outputDirFull,'xRunOutput.mat'), "stimInfo", "reconInfo");

%% Variant Summaries
%
% Build the first layer of summary figures such that you can assess the
% results of a given row of the montage (looking at the change in recon
% across constant proportion).

% Front end stuff to get proper counting
stimIndex = mod(pp, numStim);
stimIndex(stimIndex == 0) = 9;

% Grab foward cone mosaic and render matrix NOTE: Replace this with the
% loading of the stored mosaic and info without the render struct, ideally
% should be much faster. 
if (~exist(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'file'))
    error('Forward render strucure not cached')
else
    clear forwardRenderStructure;
    load(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'renderStructure');
    forwardRenderStructure = renderStructure; clear renderStructure;
    grabRenderStruct(forwardRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
        pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters);
    forwardConeMosaic = forwardRenderStructure.theConeMosaic;
end

% Where tingle is the loaded cell and 'Find me' is the name of the desired
% render output
find(strcmp(tingle, 'Find me'))

% Create a new figure if at the start of a sequence, otherwise grab the
% existing one. Then place the mosaic on it.
if stimIndex == 1
    figVariant = figure;
end
axesVariant = subplot(3,numStim,stimIndex);
forwardConeMosaic.visualizeMosaic(figVariant,axesVariant);
set(gca, 'xticklabel', [], 'yticklabel', []);
set(gca, 'xlabel', [], 'ylabel', []);

% Establish boundaries based on stimulus position
xBounds(1,:) = [pr.eccXDegs pr.eccXDegs];
yBounds(1,:) = [pr.eccYDegs pr.eccYDegs];
centerWidth = pr.stimSizeDegs/2;
xBounds(2,:) = [xBounds(1)-centerWidth, xBounds(2)+centerWidth];
yBounds(2,:) = [yBounds(1)-centerWidth, yBounds(2)+centerWidth];

% Superimpose an rectangle corrsponding to stimulus location onto the
% mosaic.
hold on
rectangle('Position', [xBounds(2,1) yBounds(2,1) (xBounds(2,2) - xBounds(2,1)) (yBounds(2,2) - yBounds(2,1))], 'LineWidth', 3)

axesVariant = subplot(3,numStim,(numStim * 2));
imshow(stimInfo.imageRGBAcrossDisplays);

axesVariant = subplot(3,numStim,(numStim * 3));
imshow(reconInfo.imageRGBAcrossDisplays);

% Double check this to see if it's sloppy code
if stimIndex == 9
    variantNames = cnv.outputDirThird(end-8:end);
    saveas(gcf, fullfile(pr.aoReconDir, pr.versEditor, cnv.generalConditions,...
        cnv.outputDirFirst, cnv.outputDirSecond, 'summaryFigs', ...
        ['SummaryRow' variantNames '.tiff']), 'tiff');
end

%% Region Proportion Summaries

%% Stim Size Summaries

end
