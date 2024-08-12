function imageInfo = buildSummaryFigs(pr, cnv, numStim, numProp, ...
    fullReconSummary, stimSummary, varargin)
% Grab image information for further processing
%
% Description:
%     Implement after running reconstruction simulations to rapidly pull
%     image results for further proccessing. Also provides the option to
%     create a low level summary figure showing the mosaics utilized, the
%     stimulus, and the reconstruction for a series of colors under one
%     proportionality.
%
% See also: aoStimReconRunMany_small_quads_chr.m, buildSummaryFigs.m
%
% History:
%   04/20/24  chr  Make it a callable function

p = inputParser;
p.addParameter('figReconRows', false, @islogical);
p.addParameter('wls', (400:1:700)', @isnumeric);
p.addParameter('allMontages', true, @islogical);
p.addParameter('scaleToMax', false, @islogical);
p.addParameter('zoomToStim', false, @islogical);
p.addParameter('centermostEW', false, @islogical);
p.addParameter('wavelengthUY', 580, @isnumeric);
p.addParameter('plotMontages', true, @islogical);
p.addParameter('plotStimvRecon', true, @islogical);
p.addParameter('plotShiftUY', true, @islogical);
p.addParameter('plotPropvRecon', true, @islogical);

parse(p, varargin{:});

close all;

%% Retrieve Recon Info
%
% Initiate directory names to descend levels for
% plotting.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
generalDir = fullfile(pr.aoReconDir, pr.versEditor, cnv.generalConditions, ...
    cnv.outputDirFirst);
generalSubDir = dir(generalDir);
imageInfo = cell(1,numStim);

% Load in the pertinent cone mosaic, note that by this procedure we
% only utilize the forward render structure information under the
% assumption that these two are similar in all meaningful ways. Need to
% reexamine this routine if that assumption fails.
if (~exist(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'file'))
    error('Forward render strucure not cached')
else
    clear forwardRenderStructure;
    load(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'renderStructure');
    forwardRenderStructure = renderStructure; clear renderStructure;
    grabRenderStruct(forwardRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
        pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters);
    theConeMosaic = forwardRenderStructure.theConeMosaic;
end

% Temporary patch to account for a previous mistake in how the cone
% proportion information was being stored. Doubles back to save the correct
% information, replace it in the render structure, and save a jpeg for
% ready access.
%
% Placing the patch here is insufficient because the whole point of loading
% out here is to do it once for one of the proportions and operate under
% the assumptions that the core underlying optics (which are all we care
% about) are held constant across the other mosaics. So essentially this
% patch only hits one of the render matrices. 
updatedMosaicConeInfo = propPatch(forwardRenderStructure.mosaicConeInfo);
forwardRenderStructure.mosaicConeInfo = updatedMosaicConeInfo;
renderStructure = forwardRenderStructure;
save(fullfile(cnv.forwardRenderDirFull, cnv.renderName),'renderStructure','-v7.3')

% Specify variables depending on the viewing display
switch (pr.viewingDisplayName)
    case 'conventional'
        viewingDisplayFieldName = 'CRT12BitDisplay';
        viewingOverwriteDisplayGamma = false;
    case 'mono'
        viewingDisplayFieldName = 'monoDisplay';
        viewingOverwriteDisplayGamma = true;
    otherwise
        error('Unknown forward display specified');
end

% Load in the display information
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
viewingDisplayLoad = load(fullfile(aoReconDir, ...
    'displays', [pr.viewingDisplayName 'Display.mat']));
eval(['viewingDisplay = viewingDisplayLoad.' viewingDisplayFieldName ';']);
clear viewingDisplayLoad

% Fix up gamma and wavelengths sampling
if (viewingOverwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    viewingDisplay.gamma = gammaOutput(:,[1 1 1]);
end

% Set up some more stuff
viewingDisplayScaleFactor = 1;
viewingDisplay = displaySet(viewingDisplay,'wave',p.Results.wls);
viewingDisplay = displaySet(viewingDisplay,'spd primaries', ...
    displayGet(viewingDisplay,'spd primaries')*viewingDisplayScaleFactor);

%% Make the Summary Montages for the reconstructions

if p.Results.plotMontages

    % An updated patch that allows cycling through the following code to
    % produce each of the desired montage possibilities. If allMontages false,
    % will produce the figures corresponding to the explicitly set scale and
    % zoom conditions. If those are not set, the default values are
    % unscaled and unzoomed
    if p.Results.allMontages
        allCombos = [true true false false; true false true false];
    else
        allCombos = [p.Results.scaleToMax; p.Results.zoomToStim];
    end

    % Then for each of the desired iterations
    for q = 1:length(allCombos)

        scaleToMax = allCombos(1,q);
        zoomToStim = allCombos(2,q);

        % Establish the tiledlayout figure and proper dimenensions
        theFig = figure();
        t = tiledlayout((numProp+1),numStim, 'TileSpacing','none');
        tileIndices = reshape(1:(numStim*(numProp+1)), numStim, (numProp+1))';
        tileIndexCounter = 1;

        fullSummary = [stimSummary; fullReconSummary];
        for j = 1:size(fullSummary, 2)
            for i = 1:size(fullSummary, 1)

                theAxes = nexttile(tileIndices(tileIndexCounter));

                if scaleToMax
                    scaleString = "Scaled";
                    meanLuminanceCdPerM2 = [];
                    [~, ~, imageLinear] = sceneFromFile( ...
                        fullSummary{i,j}.imageRGBAcrossDisplays, 'rgb', ...
                        meanLuminanceCdPerM2, theConeMosaic.Display);

                    % Determine the scale actor based on entries in the first
                    % row i.e.
                    if i == 1
                        scaleFactor = 1/max(imageLinear(:));
                    end

                    imageLinearScaled = imageLinear * scaleFactor;
                    imageLinearScaled(imageLinearScaled > 1) = 1;
                    imageRGBScaled = gammaCorrection( ...
                        imageLinearScaled, viewingDisplay);

                    if zoomToStim
                        zoomString = "Zoom";
                        imshow(imageRGBScaled(fullSummary{i,j}.idxXRange, ...
                            fullSummary{i,j}.idxXRange,:))
                    else
                        zoomString = "Full";
                        imshow(imageRGBScaled)
                    end

                else
                    scaleString = "Unscaled";
                    if zoomToStim
                        zoomString = "Zoom";
                        imshow(fullSummary{i,j}.imageRGBAcrossDisplays ...
                            (fullSummary{i,j}.idxXRange, ...
                            fullSummary{i,j}.idxXRange,:))
                    else
                        zoomString = "Full";
                        imshow(fullSummary{i,j}.imageRGBAcrossDisplays)
                    end

                end

                if j == 1 & i == 1
                    set(gca, 'xticklabel', [], 'yticklabel', []);
                    ylabel('Stim', 'FontSize',20)
                elseif j ==1
                    set(gca, 'xticklabel', [], 'yticklabel', []);
                    ylabel(sprintf('%0.2fL', pr.focalPropLList(i-1)), 'FontSize',15)
                end

                tileIndexCounter = tileIndexCounter + 1;
            end
        end

        % Spruce up the figure and save if that's the goal
        set(gcf, 'Position', [1023 7 653 970]);
        sgtitle({sprintf('Summary Montage %s %s %0.1fArcmin', ...
            scaleString, zoomString, (60*pr.stimSizeDegs))}, 'FontSize', 25);
        saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs, ...
            sprintf('summaryMontage%s%s.tiff', scaleString, zoomString)),'tiff')
    end
end

% Thoughts to add figure that shows the reconstruction montage for
% just the stimulus region where each pixel is given the mean EW value
% converted back into RGB. NOTE: This montage would not reflect the actual
% reconstruction procedure, which should include non-uniformities based on
% cone type as is currently. This is a loose approximation to better
% compare against the trends seen in the quanitification plots. 
% 
% Start w/ copy-paste of above portion to cycle through desired info cell,
% and for each chunk pull out the .meanEWFull stored in the struct. Would
% then need to convert that EW to RGB, maybe by a cycling procedure where
% begin with some combination of colors and adjust until minimize
% distinction between color and the given EW. 
% 
% Then take those values and put them into a matrix corresponding to the
% size of the stimulus region (technically can be any arbitrary size but
% this helps w/ consistency) and plot as if that was the original recon
% ouput. Should probably also scaleToMax here. Once again, note that this
% should NOT be taken to mean the reconstruction procedure returns uniform
% field reconstructions, this is an explictly set artifact. 


%% Make Stim vs Recon Plot over each of the proportions
%
% Initialize some holder variables, plot colors, and legends
reconSummaryMat = cell2mat(fullReconSummary);
stimSummaryMat = cell2mat(stimSummary);

plotColors = [(0:1/(numProp-1):1); ...
    1 - (0:1/(numProp-1):1); ...
    zeros(1, numProp)]';
plotColorsScaled = plotColors ./ max(plotColors, [], 2);

% Create the figure 2 legends based on input
legend2End = repmat(" %L", 1, numProp);
legend2Str = append(string(pr.focalPropLList*100), legend2End);
fig2Legend = num2cell(legend2Str);

if p.Results.plotStimvRecon
    theFig2 = figure();

    for i = 1:size(fullReconSummary, 1)

        % Include portion to say if centermost region or full stimulus
        % region. Default is full.
        plot([stimSummaryMat.meanEWFull], [reconSummaryMat(i,:).meanEWFull], '-o', ...
            'Color', plotColorsScaled(i,:), 'LineWidth', 3); hold on;
    end

    xlabel('Stim Wavelength', 'FontSize', 40);
    ylabel('Recon Wavelength', 'FontSize', 40);
    title(sprintf('Stim/Recon Comparison: %0.1f arcmin', (pr.stimSizeDegs * 60)), 'FontSize', 26)
    xlim([540 660])
    ylim([540 660])
    set(gcf, 'Position', [119   321   661   518]);
    box off
    axis square

    legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')

    % Save output as image and eps file for easier formatting in Adobe
    % Illustrator
    saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs,...
        sprintf('stimVsReconPlot.tiff')),'tiff');
    saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs,...
        sprintf('stimVsReconPlot.eps')),'epsc');
end

%% Make Shift in UY plot for each desired wavelength
%
% Baseline variables
stimForUYRecon = ones(1, numProp);
if p.Results.plotShiftUY

    for q = 1:length(p.Results.wavelengthUY)
        wavelengthUY = p.Results.wavelengthUY(q);

        theFig3 = figure();

        for i = 1:size(fullReconSummary, 1)

            % Capture stimulus and recon EW information in one matrix and sort
            % in ascending order.  Then make sure we have a one-to-one
            % mapping between recon and stim EW by averaging across cases
            % where there are multiple matching recon EW values.
            imageEW = double([stimSummaryMat(1,:).meanEWFull; reconSummaryMat(i,:).meanEWFull]);
            imageEWSorted = sortrows(imageEW', 2)';
            uniqueEWRecon = unique(imageEWSorted(2,:));
            uniqueEWStim = NaN*ones(size(uniqueEWRecon));
            for uu = 1:length(uniqueEWRecon)
                indexTemp = imageEWSorted(2,:) == uniqueEWRecon(uu);
                uniqueEWStim(uu) = mean(imageEWSorted(1,indexTemp));
            end

            % Perform the actual interpolation
            stimForUYRecon(i) = interp1(uniqueEWRecon,uniqueEWStim,wavelengthUY,'linear');
        end

        % Logistics check, to keep interpolation within bounds we will remove
        % any pairs that extend beyond possible stim EW range (540-680)
        plotValsUY = [pr.focalPropLList;stimForUYRecon];
        outOfRange = find(plotValsUY(2,:) < 540 | plotValsUY(2,:) > 680);
        plotValsUY(:,outOfRange) = [];

        plot(plotValsUY(1,:), plotValsUY(2,:), '-o', 'Color', ...
            plotColorsScaled(i,:), 'Linewidth', 5);

        xlabel('Proportion L', 'FontSize', 40);
        ylabel('Stim Wavelength', 'FontSize', 40);
        title(sprintf('Shift in UY: %d nm %0.1f arcmin', wavelengthUY,  ...
            (pr.stimSizeDegs * 60)), 'FontSize', 26)
        xlim([0 1])
        ylim([540 660])
        set(gcf, 'Position', [119   321   661   518]);
        box off
        axis square

        % Save output as image and eps file for easier formatting in Adobe
        % Illustrator
        saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs,...
            sprintf('shiftUYPlot%d.tiff', wavelengthUY)),'tiff');
        saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs,...
            sprintf('shiftUYPlot%d.eps', wavelengthUY)),'epsc');
    end
end

%% Make Prop vs Recon Plot over each of the stimuli
%
% Initialize a new set of colors
plotColors = [(0:1/(numStim-1):1); ...
    1 - (0:1/(numStim-1):1); ...
    zeros(1, numStim)]';
plotColorsScaled = plotColors ./ max(plotColors, [], 2);

% Create the figure 3 legends based on input
legend3End = repmat(" nm", 1, numStim);
legend3Str = append(string([stimSummaryMat(1,:).meanEWFull]), legend3End);
fig3Legend = num2cell(legend3Str);

if p.Results.plotPropvRecon
    theFig4 = figure();

    for i = 1:size(fullReconSummary, 2)
        % Include portion to say if centermost region or full stimulus
        % region. Default is full.
        plot(pr.focalPropLList, [reconSummaryMat(:,i).meanEWFull], '-o', ...
            'Color', plotColorsScaled(i,:), 'LineWidth', 3); hold on;
    end

    xlabel('Proportion L', 'FontSize', 40);
    ylabel('Recon Wavelength', 'FontSize', 40);
    title(sprintf('Prop/Recon Comparison: %0.1f arcmin', (pr.stimSizeDegs * 60)), 'FontSize', 26)
    xlim([0 1])
    ylim([540 680])
    set(gcf, 'Position', [119   321   661   518]);
    box off
    axis square

    legend(fig3Legend, 'NumColumns',2, 'Location', 'northwest')

    % Save output as image and eps file for easier formatting in Adobe
    % Illustrator
    saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs,...
        sprintf('propVsReconPlot.tiff')),'tiff');
    saveas(gcf, fullfile(cnv.outputSubdirSummaryFigs,...
        sprintf('propVsReconPlot.eps')),'epsc');
end
end
