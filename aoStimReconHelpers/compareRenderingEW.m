function imageEW = compareRenderingEW(stimImageRGBFormer, reconImageRGBFormer, ...
stimImagergbLinear, reconImagergbLinear, startDisplayName, viewingDisplayName, ...
idxXRange, varargin)

% This is originally the ReconRerunFigs program. The idea is that it
% needs to be cleaned up and restructured, and then fold back into what
% has become the RenderRecons document.
%
% Description:
%    Collect up a set of reconstructions and render/analyze.
%
%    Patterned after aoStimReconRerunFigs, but simplifying.
%
% See also: aoStimRecon, aoStimReconRunMany, correctForViewing

% Notes:
%    Cone ratios are such that the first number in a series (i.e. 110) is 90% L and then
%    systematically goes down towards 10%L for a given chunk.
%
%    Also I can double check as for the >136 I may not have pushed those edits to the git. 137-140 are positional changes
%    in same proportionality, 141 is all L surround, 142 is all M surround

% History:
%   06/01/23  chr  Organize into its own file
%   10/19/23  dhb  Wrote from chr program.
%   01/21/24  chr  Set up this script for processing of wavelength values

%Random 


% Parse key value pairs
p = inputParser;
p.addParameter('showFigs', false, @islogical);
p.addParameter('inwardMove', 0, @isnumeric)
p.addParameter('linearInput', false, @islogical);
p.addParameter('viewingDisplayScaleFactor',1,@isnumeric);
p.addParameter('wls',(400:10:700)',@isnumeric);
p.addParameter('verbose',false,@islogical);
p.addParameter('SRGB',false,@islogical);
p.addParameter('scaleToMax',false,@islogical)
parse(p, varargin{:});

% Establish a struct for the EW output values to then be analyzed
imageEW = struct;

% Pull the stimulus and recon image used in the montages, note the stimulus
% has already been corrected from mono to conventional dislay based on the
% approach previously taken by the code
% stimImageRGBFormer = cfvStim.stimulusRGBScaled{1};
% reconImageRGBFormer = cfvRecon.reconScaledRGB{1};

% Want to overhaul since still feel like missing the point at some instance
% prior
% inputImageRGB = stimulusImageRGB;
% startDisplayName = rrf.startDisplayName;
% viewingDisplayName = rrf.viewingDisplayName;
% viewingDisplayScaleFactor = rrf.stimDispScale;
% SRGB = false;


% Specify variables depending on the start display
switch (startDisplayName)
    case 'conventional'
        startDisplayFieldName = 'CRT12BitDisplay';
        startOverwriteDisplayGamma = false;
    case 'mono'
        startDisplayFieldName = 'monoDisplay';
        startOverwriteDisplayGamma = true;
    otherwise
        error('Unknown recon display specified');
end

% Specify variables depending on the viewing display
switch (viewingDisplayName)
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
startDisplayLoad = load(fullfile(aoReconDir, 'displays', [startDisplayName 'Display.mat']));
eval(['startDisplay = startDisplayLoad.' startDisplayFieldName ';']);
viewingDisplayLoad = load(fullfile(aoReconDir, 'displays', [viewingDisplayName 'Display.mat']));
eval(['viewingDisplay = viewingDisplayLoad.' viewingDisplayFieldName ';']);
clear startDisplayLoad viewingDisplayLoad


% Sanity check by redoing the display correction based on the updates,
% utilizes the RenderAcrossDisplay functionality and ideally should be the
% same as before
stimImageRGBUncorrected = gammaCorrection(stimImagergbLinear,startDisplay);
[stimImageRGBCurrent,stimImagergbTruncated,stimImagergb] = RGBRenderAcrossDisplays(stimImageRGBUncorrected, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB);

reconImageRGBUncorrected = gammaCorrection(reconImagergbLinear,startDisplay);
[reconImageRGBCurrent,reconImagergbTruncated,reconImagergb] = RGBRenderAcrossDisplays(reconImageRGBUncorrected, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB);



% Adjust bounds inward in an attempt to avoid skewing by edge artifacts
% seen in reconstructions. Could pose issues if image is too small. 
idxXRangeNew = (idxXRange(1) + p.Results.inwardMove) : (idxXRange(end) - p.Results.inwardMove);



% Calculate equivalent wavelength for the forward/recon pair and the
% mono/conventional pair
stimImageRGBforEW1 = (stimImageRGBFormer(idxXRangeNew, idxXRangeNew, :));
reconImageRGBforEW1 = (reconImageRGBFormer(idxXRangeNew, idxXRangeNew, :));

stimImageRGBforEW2 = (stimImageRGBCurrent(idxXRangeNew, idxXRangeNew, :));
reconImageRGBforEW2 = (reconImageRGBCurrent(idxXRangeNew, idxXRangeNew, :));

stimImageRGBforEWUncorrected = (stimImageRGBUncorrected(idxXRangeNew, idxXRangeNew, :));
reconImageRGBforEWUncorrected = (reconImageRGBUncorrected(idxXRangeNew, idxXRangeNew, :));


% Calculate the EW on the stim/recon pair formerly corrected to viewing
% display
[imageEW.stimEWFormer] = ...
    RGBToEquivWavelength(stimImageRGBforEW1, viewingDisplay);
[imageEW.reconEWFormer] = ...
    RGBToEquivWavelength(reconImageRGBforEW1, viewingDisplay);

% Calculate the EW on the stim/recon pair currently corrected to viewing
% display
[imageEW.stimEWCurrent] = ...
    RGBToEquivWavelength(stimImageRGBforEW2, viewingDisplay);
[imageEW.reconEWCurrent] = ...
    RGBToEquivWavelength(reconImageRGBforEW2, viewingDisplay);

% Calculate the EW on the stim/recon pair not corrected to viewing
% display
[imageEW.stimEWUncorrected] = ...
    RGBToEquivWavelength(stimImageRGBforEWUncorrected, startDisplay);
[imageEW.reconEWUncorrected] = ...
    RGBToEquivWavelength(reconImageRGBforEWUncorrected, startDisplay);



if p.Results.showFigs
    figure()
    subplot(3,2,1)
    imshow(stimImageRGBforEW1)
    title(sprintf('Mono Stim EW Corrected 1: %d', ...
        int64(mean(imageEW.stimEWFormer, 'all'))));

    subplot(3,2,2)
    imshow(reconImageRGBforEW1)
    title(sprintf('Mono Recon EW Corrected 1: %d', ...
        int64(mean(imageEW.reconEWFormer, 'all'))));

    subplot(3,2,3)
    imshow(stimImageRGBforEW2)
    title(sprintf('Mono Stim EW Corrected 2: %d', ...
        int64(mean(imageEW.stimEWCurrent, 'all'))));

    subplot(3,2,4)
    imshow(reconImageRGBforEW2)
    title(sprintf('Mono Recon EW Corrected 2: %d', ...
        int64(mean(imageEW.reconEWCurrent, 'all'))));

    subplot(3,2,5)
    imshow(stimImageRGBforEWUncorrected)
    title(sprintf('Conv Stim EW Uncorrected 2: %d', ...
        int64(mean(imageEW.stimEWUncorrected, 'all'))));

    subplot(3,2,6)
    imshow(reconImageRGBforEWUncorrected)
    title(sprintf('Conv Recon EW Uncorrected 2: %d', ...
        int64(mean(imageEW.reconEWUncorrected, 'all'))));
end



%
%
% % Pull the stimulus and recon image used in the montages, note the stimulus
% % has already been corrected from mono to conventional display
% inputImageRGB = cfvStim.stimulusRGBScaled{1};
% reconImageRGB = cfvRecon.reconScaledRGB{1};
%
% % Want to overhaul since still feel like missing the point at some instance
% % prior
% % inputImageRGB = stimulusImageRGB;
% startDisplayName = rrf.startDisplayName;
% viewingDisplayName = rrf.viewingDisplayName;
% viewingDisplayScaleFactor = rrf.stimDispScale;
% SRGB = false;
%
% % Specify variables depending on the start display
% switch (startDisplayName)
%     case 'conventional'
%         startDisplayFieldName = 'CRT12BitDisplay';
%         startOverwriteDisplayGamma = false;
%     case 'mono'
%         startDisplayFieldName = 'monoDisplay';
%         startOverwriteDisplayGamma = true;
%     otherwise
%         error('Unknown recon display specified');
% end
%
% % Specify variables depending on the viewing display
% switch (viewingDisplayName)
%     case 'conventional'
%         viewingDisplayFieldName = 'CRT12BitDisplay';
%         viewingOverwriteDisplayGamma = false;
%     case 'mono'
%         viewingDisplayFieldName = 'monoDisplay';
%         viewingOverwriteDisplayGamma = true;
%     otherwise
%         error('Unknown forward display specified');
% end
%
%
% % Load in the display information
% aoReconDir = getpref('ISETImagePipeline','aoReconDir');
% viewingDisplayLoad = load(fullfile(aoReconDir, 'displays', [viewingDisplayName 'Display.mat']));
% eval(['viewingDisplay = viewingDisplayLoad.' viewingDisplayFieldName ';']);
% clear startDisplayLoad viewingDisplayLoad
%
% idxLBN = idxLB + 3;
% idxUBN = idxUB - 3;
%
% % Calculate equivalent wavelength for the forward/recon pair and the
% % mono/conventional pair
% inputImageRGBforEW = (inputImageRGB(idxLBN:idxUBN, idxLBN:idxUBN, :));
% reconImageRGBforEW = (reconImageRGB(idxLBN:idxUBN, idxLBN:idxUBN, :));
%
% % Calculate the EW on the stim/recon pair already corrected to viewing
% % display
% [equivWavelengthInput] = ...
%     RGBToEquivWavelength(inputImageRGBforEW, viewingDisplay);
% [equivWavelengthRecon] = ...
%     RGBToEquivWavelength(reconImageRGBforEW, viewingDisplay);
