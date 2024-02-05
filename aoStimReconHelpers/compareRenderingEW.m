function imageEW = compareRenderingEW(stimImageRGBFormer, reconImageRGBFormer, ...
stimImagergbLinear, reconImagergbLinear, startDisplayName, viewingDisplayName, ...
idxXRange, varargin)

% This is originally the ReconRerunFigs program. The idea is that it
% needs to be cleaned up and restructured, and then fold back into what
% has become the RenderRecons document.
%
% Description:
%    Fed into aoStimReconRerunFigs for comparison across EW calculations
%
% See also: aoStimRecon, aoStimReconRunMany, correctForViewing

% Notes:
%    Cone ratios are such that the first number in a series (i.e. 110) is 90% L and then
%    systematically goes down towards 10%L for a given chunk.
%
%    Also I can double check as for the >136 I may not have pushed those edits to the git. 137-140 are positional changes
%    in same proportionality, 141 is all L surround, 142 is all M surround

% History:
%   01/21/24  chr  Set up this script for processing of wavelength values

% Parse key value pairs
p = inputParser;
p.addParameter('showFigs', false, @islogical);
p.addParameter('inwardMove', 0, @isnumeric)
p.addParameter('linearInput', false, @islogical);
p.addParameter('viewingDisplayScaleFactor',1,@isnumeric);
p.addParameter('wls',(400:1:700)',@isnumeric);
p.addParameter('verbose',false,@islogical);
p.addParameter('SRGB',false,@islogical);
p.addParameter('scaleToMax',false,@islogical)
parse(p, varargin{:});

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


% Establish a struct for the EW output values to then be analyzed
imageEW = struct;
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
displayGammaBits = 12;
displayGammaGamma = 2;
% overwriteDisplayGamma = true;
viewingDisplayScaleFactor = 1;


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
% aoReconDir = getpref('ISETImagePipeline','aoReconDir');
startDisplayLoad = load(fullfile(aoReconDir, 'displays', [startDisplayName 'Display.mat']));
eval(['startDisplay = startDisplayLoad.' startDisplayFieldName ';']);
viewingDisplayLoad = load(fullfile(aoReconDir, 'displays', [viewingDisplayName 'Display.mat']));
eval(['viewingDisplay = viewingDisplayLoad.' viewingDisplayFieldName ';']);
clear startDisplayLoad viewingDisplayLoad

% Fix up gamma and wavelengths sampling
if (startOverwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    startDisplay.gamma = gammaOutput(:,[1 1 1]);
end

if (viewingOverwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    viewingDisplay.gamma = gammaOutput(:,[1 1 1]);
end
    


% wls = (400:1:700)';
startDisplay = displaySet(startDisplay,'wave',p.Results.wls);
viewingDisplay = displaySet(viewingDisplay,'wave',p.Results.wls);
viewingDisplay = displaySet(viewingDisplay,'spd primaries',displayGet(viewingDisplay,'spd primaries')*viewingDisplayScaleFactor);


% Sanity check by redoing the display correction based on the updates,
% utilizes the RenderAcrossDisplay functionality and ideally should be the
% same as before
stimImageRGBUncorrected = gammaCorrection(stimImagergbLinear,startDisplay);
[stimImageRGBCurrent,stimImagergbTruncated,stimImagergb] = RGBRenderAcrossDisplays(stimImageRGBUncorrected, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB, ...
    'wls', p.Results.wls);

reconImageRGBUncorrected = gammaCorrection(reconImagergbLinear,startDisplay);
[reconImageRGBCurrent,reconImagergbTruncated,reconImagergb] = RGBRenderAcrossDisplays(reconImageRGBUncorrected, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB, ...
    'wls', p.Results.wls);



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
    imshow(stimImageRGBFormer)
    title(sprintf('Mono Stim EW Corrected 1: %d', ...
        int64(mean(imageEW.stimEWFormer, 'all'))));

    subplot(3,2,2)
    imshow(reconImageRGBFormer)
    title(sprintf('Mono Recon EW Corrected 1: %d', ...
        int64(mean(imageEW.reconEWFormer, 'all'))));

    subplot(3,2,3)
    imshow(stimImageRGBCurrent)
    title(sprintf('RenderAcrossDisplays Conv Stim: %d', ...
        int64(mean(imageEW.stimEWCurrent, 'all'))));

    subplot(3,2,4)
    imshow(reconImageRGBCurrent)
    title(sprintf('RenderAcrossDisplays Conv Recon: %d', ...
        int64(mean(imageEW.reconEWCurrent, 'all'))));

    subplot(3,2,5)
    imshow(stimImageRGBUncorrected)
    title(sprintf('Mono Stim: %d', ...
        int64(mean(imageEW.stimEWUncorrected, 'all'))));

    subplot(3,2,6)
    imshow(reconImageRGBUncorrected)
    title(sprintf('Mono Recon: %d', ...
        int64(mean(imageEW.reconEWUncorrected, 'all'))));
end



imageEW.stimImageRGB = stimImageRGBCurrent; % Post Disp Correction
imageEW.reconImageRGB = reconImageRGBCurrent; % Post Disp Correction



% Scale to max if specified
if (p.Results.scaleToMax)
    stimImagergbLinearScaled = stimImagergbLinear / max(stimImagergbLinear(:));
    stimImageRGBUncorrected = gammaCorrection(stimImagergbLinearScaled,startDisplay);

    reconImagergbLinearScaled = reconImagergbLinear / max(reconImagergbLinear(:));
    reconImageRGBUncorrected = gammaCorrection(reconImagergbLinearScaled,startDisplay);
end

imageEW.stimImageRGBTest = stimImageRGBUncorrected; % Prior to Disp Correction
imageEW.reconImageRGBTest = reconImageRGBUncorrected; % Prior to Disp Correction