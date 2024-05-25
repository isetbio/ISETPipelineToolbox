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
p.addParameter('scaleToMax', false, @islogical);
p.addParameter('wls', (400:1:700)', @isnumeric);
p.addParameter('zoomToStim', false, @islogical);
p.addParameter('centermostEW', false, @islogical);
parse(p, varargin{:});

close all;

%% Retrieve Recon Info
%
% Initiate directory names to descend levels for plotting.
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

%% Make the recon image montage

testing = true;

if testing
    
    % Establish the tiledlayout figure and proper dimenensions
    theFig = figure();
    t = tiledlayout((numProp+1),numStim, 'TileSpacing','none');
    tileIndices = reshape(1:(numStim*(numProp+1)), numStim, (numProp+1))';
    tileIndexCounter = 1;
    
    fullSummary = [stimSummary; fullReconSummary];
    for j = 1:size(fullSummary, 2)
        for i = 1:size(fullSummary, 1)
            
            theAxes = nexttile(tileIndices(tileIndexCounter));
            
            if p.Results.scaleToMax
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
                
                if p.Results.zoomToStim
                    zoomString = "Zoom";
                    imshow(imageRGBScaled(fullSummary{i,j}.idxXRange, ...
                        fullSummary{i,j}.idxXRange,:))
                else
                    zoomString = "Full";
                    imshow(imageRGBScaled)
                end
                
            else
                scaleString = "Unscaled";
                if p.Results.zoomToStim
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
    saveas(gcf, fullfile(generalDir, cnv.outputDirSecond, ...
        sprintf('summaryMontage%s%s.tiff', scaleString, zoomString)),'tiff')
end

%% Make the Summary plots

john = true;

if john
    plotColors = [(0:1/(numProp-1):1); ...
        1 - (0:1/(numProp-1):1); ...
        zeros(1, numProp)]';
    plotColorsScaled = plotColors ./ max(plotColors, [], 2);
    
    % Create the figure 2 legends based on input
    legend2End = repmat(" %L", 1, numProp);
    legend2Str = append(string(pr.focalPropLList*100), legend2End);
    fig2Legend = num2cell(legend2Str);
    
    berryRecon = cell2mat(fullReconSummary);
    berryStim = cell2mat(stimSummary);
    theFig2 = figure();
    
    for i = 1:size(fullReconSummary, 1)
        
        % Include portion to say if centermost region or full stimulus
        % region. Default is full. 
        plot([berryStim.meanEWFull], [berryRecon(i,:).meanEWFull], '-o', ...
            'Color', plotColorsScaled(i,:), 'LineWidth', 3); hold on;
    end
    
    xlabel('Stim Wavelength', 'FontSize', 40);
    ylabel('Recon Wavelength', 'FontSize', 40);
    title(sprintf('Stim/Recon Comparison: %d arcmin', (pr.stimSizeDegs * 60)), 'FontSize', 26)
    xlim([540 660])
    ylim([540 660])
    set(gcf, 'Position', [119   321   661   518]);
    box off
    axis square
    
    legend(fig2Legend, 'NumColumns',2, 'Location', 'northwest')
    
    % Save output as image and eps file for easier formatting in Adobe
    % Illustrator
    saveas(gcf, fullfile(generalDir, cnv.outputDirSecond,...
        sprintf('stimVsReconPlot.tiff')),'tiff');
    
    saveas(gcf, fullfile(generalDir, cnv.outputDirSecond,...
        sprintf('stimVsReconPlot.eps')),'epsc');
    
end

jerry = true;

if jerry
end


end
