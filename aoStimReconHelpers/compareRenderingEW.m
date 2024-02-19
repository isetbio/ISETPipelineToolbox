function imageEW = compareRenderingEW(stimImageRGBFormer, reconImageRGBFormer, ...
stimImagergbLinear, reconImagergbLinear, startDisplayName, viewingDisplayName, ...
idxXRange, varargin)
% Facilitate image display correction and Equivalent Wavelength calculation
%
% imageEW = compareRenderingEW[stimImageRGBFormer, reconImageRGBFormer, 
% stimImagergbLinear, reconImagergbLinear, startDisplayName, viewingDisplayName,
% idxXRange, varagin]
%
% Description:
%    This is an organizational wrapper function to facilitate the rendering of
%    stimulus/reconstruction pairs across different displays (i.e. mono to
%    conventional). The first two arguments are image inputs from previous
%    iterations incorporated in here as a sanity check. However, going
%    forward will be conducting the correction here post-simulation
%    instead. 
% 
%    Also facilitates the calculation of equivalent wavelength for
%    stim/recon pairs. Should emphasize that this calculation should be
%    conducted using the numerical values prior to any display correction
%    (i.e. if correcting from mono to conventional, use the mono values to 
%    get EW). Various other methods are used to calculate EW in the script
%    below but this is purely to regain intuition about values, the only
%    one of actual concern is that saved to the final imageEW output struct
% 
%    This function is called within the larger wrapper script
%    aoStimReconRerunFigs_chr.m which applies applies this to an entire
%    directory of simulation output files. 
%
%
% See also: RGBRenderAcrossDisplays, RGBToEquivalentWavelength

% Notes:
%    Cone ratios are such that the first number in a series (i.e. 110) is 90% L and then
%    systematically goes down towards 10%L for a given chunk.
%
%    Also I can double check as for the >136 I may not have pushed those edits to the git. 137-140 are positional changes
%    in same proportionality, 141 is all L surround, 142 is all M surround

% History:
%   01/21/24  chr  Set up this script for processing of wavelength values
%   02/19/24  chr  Cleaning

% Parse key value pairs
p = inputParser;
p.addParameter('showFigs', false, @islogical);
p.addParameter('inwardMove', true, @islogical)
p.addParameter('linearInput', false, @islogical);
p.addParameter('viewingDisplayScaleFactor',1,@isnumeric);
p.addParameter('wls',(400:1:700)',@isnumeric);
p.addParameter('verbose',false,@islogical);
p.addParameter('SRGB',false,@islogical);
p.addParameter('scaleToMax',false,@islogical)
parse(p, varargin{:});

% Establish a struct for the EW output values to then be analyzed
imageEW = struct;
aoReconDir = getpref('ISETImagePipeline','aoReconDir');
displayGammaBits = 12;
displayGammaGamma = 2;
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
stimImageRGBGamma = gammaCorrection(stimImagergbLinear,startDisplay);
[stimImageRGBGammaDisplay,stimImagergbDisplayTruncated,stimImagergbDisplay] = ...
    RGBRenderAcrossDisplays(stimImageRGBGamma, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB, ...
    'wls', p.Results.wls);

reconImageRGBGamma = gammaCorrection(reconImagergbLinear,startDisplay);
[reconImageRGBGammaDisplay,reconImagergbDisplayTruncated,reconImagergbDisplay] = ...
    RGBRenderAcrossDisplays(reconImageRGBGamma, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB, ...
    'wls', p.Results.wls);

% If inwardMove set to true, establish a central patch for calculation
% of the Equivalent Wavelength to avoid edge artifact impacting data
if p.Results.inwardMove 
    if isodd(length(idxXRange))
        centerPoint = ceil(length(idxXRange)/2);
        centerSpread = floor(length(idxXRange) / 4);
        idxXRangeNew = idxXRange(centerPoint - centerSpread: centerPoint + centerSpread);
    else
        centerPoint = ceil(length(idxXRange)/2);
        centerSpread = floor(length(idxXRange) / 4);stimEWFormer
        idxXRangeNew = idxXRange(centerPoint - centerSpread + 1: centerPoint + centerSpread);
    end
end

% Apply the new central patch to the stimulus and reconstruction
stimImageRGBforEWOldVersion = (stimImageRGBFormer(idxXRangeNew, idxXRangeNew, :));
reconImageRGBforEWOldVersion = (reconImageRGBFormer(idxXRangeNew, idxXRangeNew, :));

stimImageRGBforEWDisplay = (stimImageRGBGammaDisplay(idxXRangeNew, idxXRangeNew, :));
reconImageRGBforEWDisplay = (reconImageRGBGammaDisplay(idxXRangeNew, idxXRangeNew, :));

stimImageRGBforEW = (stimImageRGBGamma(idxXRangeNew, idxXRangeNew, :));
reconImageRGBforEW = (reconImageRGBGamma(idxXRangeNew, idxXRangeNew, :));

% Calculate the EW for the old version of stim/recon, after display 
% correction but prior to overhaul
[stimEWOldVersion] = ...
    RGBToEquivWavelength(stimImageRGBforEWOldVersion, viewingDisplay);
[reconEWOldVersion] = ...
    RGBToEquivWavelength(reconImageRGBforEWOldVersion, viewingDisplay);

% Calculate the EW for the current version of stim/recon, after display
% correction
[stimEWDisplay] = ...
    RGBToEquivWavelength(stimImageRGBforEWDisplay, viewingDisplay);
[reconEWDisplay] = ...
    RGBToEquivWavelength(reconImageRGBforEWDisplay, viewingDisplay);

% Calculate the EW for current version of stim/recon, prior to display
% correction. NOTE! This is the value that should be used for reporting. EW
% should be calculated prior to correction, images should be viewed after
[imageEW.stimEW] = ...
    RGBToEquivWavelength(stimImageRGBforEW, startDisplay);
[imageEW.reconEW] = ...
    RGBToEquivWavelength(reconImageRGBforEW, startDisplay);

% Summary figure as a sanity check just to regain intuition about certain
% aspects of the calculations. 
if p.Results.showFigs
    figure()
    subplot(3,2,1)
    imshow(stimImageRGBFormer)
    title(sprintf('Mono Stim EW Corrected 1: %d', ...
        int64(mean(stimEWOldVersion, 'all'))));

    subplot(3,2,2)
    imshow(reconImageRGBFormer)
    title(sprintf('Mono Recon EW Corrected 1: %d', ...
        int64(mean(reconEWOldVersion, 'all'))));

    subplot(3,2,3)
    imshow(stimImageRGBGammaDisplay)
    title(sprintf('RenderAcrossDisplays Conv Stim: %d', ...
        int64(mean(stimEWDisplay, 'all'))));

    subplot(3,2,4)
    imshow(reconImageRGBGammaDisplay)
    title(sprintf('RenderAcrossDisplays Conv Recon: %d', ...
        int64(mean(reconEWDisplay, 'all'))));

    subplot(3,2,5)
    imshow(stimImageRGBGamma)
    title(sprintf('Mono Stim: %d', ...
        int64(mean(imageEW.stimEW, 'all'))));

    subplot(3,2,6)
    imshow(reconImageRGBGamma)
    title(sprintf('Mono Recon: %d', ...
        int64(mean(imageEW.reconEW, 'all'))));
end

% Output a version of the images after rendering across Displays
imageEW.stimImageRGB = stimImageRGBGammaDisplay; 
imageEW.reconImageRGB = reconImageRGBGammaDisplay; 