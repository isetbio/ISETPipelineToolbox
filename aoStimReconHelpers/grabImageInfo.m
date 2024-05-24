function [stimSummary, reconSummary] = grabImageInfo(pr, cnv, numStim, varargin)
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
parse(p, varargin{:});

close all;

%% Retrieve Recon Info
%
% Initiate directory names to descend levels for plotting.
generalDir = fullfile(pr.aoReconDir, pr.versEditor, cnv.generalConditions, ...
    cnv.outputDirFirst, cnv.outputDirSecond, cnv.outputDirThird);
generalSubDir = dir(generalDir);
stimSummary = cell(1,numStim);
reconSummary = cell(1,numStim);
infoCounter = 1;

% Account for some of the subdirectories created by default to avoid errors
% in starting position when going down the loop. Another approach to
% consider would be to ID the subdirectories with names that match the
% desired format (i.e. stimColor) but this seems to do the job.
if isempty(generalSubDir)
    error(['No output information stored at desired directory: ', generalDir])
end

subDirNames = extractfield(generalSubDir, 'name');
subDirNamesFullString = cell2mat(subDirNames);
if contains(subDirNamesFullString, '.DS_Store')
    startDirIndex = 4;
else
    startDirIndex = 3;
end

% If we don't want to visualize/update the recon plots then simply pull the
% pertinent image info as output
if ~p.Results.figReconRows
    
    % For each of the contained subdirectories
    for i = startDirIndex:length(generalSubDir)
        if generalSubDir(i).isdir
            
            % Load the ouptut file to get the stim and recon images, store
            % this information for future use
            load(fullfile(generalDir, generalSubDir(i).name, ...
                'xRunOutput.mat'), "stimInfo", "reconInfo", ...
                "idxXRange");
            
            % Include information on just the stim location for figures.
            % Operates under the assumption that stimuli are squares.
            stimInfo.idxXRange = idxXRange;
            reconInfo.idxXRange = idxXRange;
            
            % Place the pertinent info into a cell for ease of access
            stimSummary(infoCounter) = {stimInfo};
            reconSummary(infoCounter) = {reconInfo};
            infoCounter = infoCounter + 1;
        end
    end
    
    % Otherwise do the same but also build the actual figure itself.
else
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
    
    % Establish the tiledlayout figure and proper dimenensions
    theFig = figure();
    t = tiledlayout(3,numStim, 'TileSpacing','tight');
    tileIndices = reshape(1:(3*numStim), numStim, 3)';
    tileIndexCounter = 1;
    
    % Set the boundaries based on stimulus position to superimpose on the
    % mosaic for viewing
    xBounds(1,:) = [pr.eccXDegs pr.eccXDegs];
    yBounds(1,:) = [pr.eccYDegs pr.eccYDegs];
    centerWidth = pr.stimSizeDegs/2;
    xBounds(2,:) = [xBounds(1)-centerWidth, xBounds(2)+centerWidth];
    yBounds(2,:) = [yBounds(1)-centerWidth, yBounds(2)+centerWidth];
    
    % For each of the contained subdirectories
    for i = startDirIndex:length(generalSubDir)
        if generalSubDir(i).isdir
            
            % Load the ouptut file to get the stim and recon images, store
            % this information for future use
            load(fullfile(generalDir, generalSubDir(i).name, ...
                'xRunOutput.mat'), "stimInfo", "reconInfo", ...
                "idxXRange");
            
            % Include information on just the stim location for figures
            % Operates under the assumption that stimuli are squares.
            stimInfo.idxXRange = idxXRange;
            reconInfo.idxXRange = idxXRange;
            
            % Place the pertinent info into a cell for ease of access
            stimSummary(infoCounter) = {stimInfo};
            reconSummary(infoCounter) = {reconInfo};
            infoCounter = infoCounter + 1;
            
            % If it's a directory, start by plotting the mosaic
            theAxes = nexttile(tileIndices(tileIndexCounter));
            theConeMosaic.visualizeMosaic(theFig,theAxes);
            set(gca, 'xticklabel', [], 'yticklabel', []);
            set(gca, 'xlabel', [], 'ylabel', []);
            [~, plotColm] = find(tileIndices == tileIndices(tileIndexCounter));
            if plotColm == 1
                ylabel('Mosaic')
            end
            tileIndexCounter = tileIndexCounter+1;
            
            % Outline the stimulus region
            hold on
            rectangle('Position', [xBounds(2,1) yBounds(2,1) ...
                (xBounds(2,2) - xBounds(2,1)) ...
                (yBounds(2,2) - yBounds(2,1))], 'LineWidth', 3)
            
            % Plot an image of the stimulus in the second row, apply
            % scaling if desired
            nexttile(tileIndices(tileIndexCounter))
            if p.Results.scaleToMax
                scaleString = "Scaled";
                meanLuminanceCdPerM2 = [];
                [~, ~, stimImageLinear] = sceneFromFile( ...
                    stimInfo.imageRGBAcrossDisplays, 'rgb', ...
                    meanLuminanceCdPerM2, theConeMosaic.Display);
                scaleFactor = 1/max(stimImageLinear(:));
                stimImageLinearScaled = stimImageLinear * scaleFactor;
                stimImageLinearScaled(stimImageLinearScaled > 1) = 1;
                stimImageRGBScaled = gammaCorrection( ...
                    stimImageLinearScaled, viewingDisplay);
                
                if p.Results.zoomToStim
                    zoomString = "Zoom";
                    imshow(stimImageRGBScaled(stimInfo.idxXRange, ...
                        stimInfo.idxXRange, :))
                else
                    zoomString = "Full";
                    imshow(stimImageRGBScaled)
                end
                
            else
                scaleString = "Unscaled";
                
                if p.Results.zoomToStim
                    zoomString = "Zoom";
                    imshow(stimInfo.imageRGBAcrossDisplays(...
                        stimInfo.idxXRange, stimInfo.idxXRange, ...
                        :))
                else
                    zoomString = "Full";
                    imshow(stimInfo.imageRGBAcrossDisplays);
                end
                
            end
            
            % Determine if image is in the first column for labeling
            % purposes
            [~, plotColm] = find(tileIndices == tileIndices(tileIndexCounter));
            if plotColm == 1
                set(gca, 'xticklabel', [], 'yticklabel', []);
                ylabel('Stimulus', 'FontSize', 20)
            end
            tileIndexCounter = tileIndexCounter+1;
            
            % Plot an image of the reconstruction in the third row, apply
            % scaling if desired
            nexttile(tileIndices(tileIndexCounter))
            if p.Results.scaleToMax
                meanLuminanceCdPerM2 = [];
                [~, ~, reconImageLinear] = sceneFromFile( ...
                    reconInfo.imageRGBAcrossDisplays, 'rgb', ...
                    meanLuminanceCdPerM2, theConeMosaic.Display);
                reconImageLinearScaled = reconImageLinear * scaleFactor;
                reconImageLinearScaled(reconImageLinearScaled > 1) = 1;
                reconImageRGBScaled = gammaCorrection( ...
                    reconImageLinearScaled, viewingDisplay);
                
                if p.Results.zoomToStim
                    imshow(reconImageRGBScaled(reconInfo.idxXRange, ...
                        reconInfo.idxXRange, :))
                else
                    imshow(reconImageRGBScaled)
                end
                
            else
                scaleString = "Unscaled";
                
                if p.Results.zoomToStim
                    imshow(reconInfo.imageRGBAcrossDisplays(...
                        reconInfo.idxXRange, reconInfo.idxXRange, ...
                        :))
                else
                    imshow(reconInfo.imageRGBAcrossDisplays);
                end
            end
            
                % Determine if image is in the first column for labeling
                % purposes
                [~, plotColm] = find(tileIndices == tileIndices(tileIndexCounter));
                if plotColm == 1
                    set(gca, 'xticklabel', [], 'yticklabel', []);
                    ylabel('Recon', 'FontSize',20)
                end
                tileIndexCounter = tileIndexCounter+1;
            end
        end
        
        % Spruce up the figure and save if that's the goal
        set(gcf, 'Position', [297 40 1601 633]);
        sgtitle({sprintf('Recons %0.2fL %s %s %0.1fArcmin', ...
            pr.focalPropL, scaleString, zoomString, (60*pr.stimSizeDegs))}, ...
            'FontSize', 30);
        saveas(gcf, fullfile(generalDir, ...
            sprintf('summaryFig%s%s.tiff', scaleString, zoomString)),'tiff')
    end
end