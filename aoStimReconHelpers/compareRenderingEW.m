function imageInfo = compareRenderingEW(imagergbLinear, startDisplayName, ...
    viewingDisplayName, idxXRange, varargin)
% Facilitate image display correction and Equivalent Wavelength calculation
%
% imageEW = compareRenderingEW[stimImageRGBFormer, reconImageRGBFormer,
% inputImagergbLinear, reconImagergbLinear, startDisplayName, viewingDisplayName,
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
%   03/29/24  chr  Adjust function so serves more as a wrapper 

% Parse key value pairs
p = inputParser;
p.addParameter('showFigs', false, @islogical);
p.addParameter('useFullEW', true, @islogical)
p.addParameter('linearInput', false, @islogical);
p.addParameter('viewingDisplayScaleFactor',1,@isnumeric);
p.addParameter('wls',(400:1:700)',@isnumeric); %%%%%
p.addParameter('verbose',false,@islogical);
p.addParameter('SRGB',false,@islogical);
p.addParameter('scaleToMax',false,@islogical)
parse(p, varargin{:});

% Establish a struct for the EW output values to then be analyzed
imageInfo = struct;
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

% Set up some more stuff
startDisplay = displaySet(startDisplay,'wave',p.Results.wls);
viewingDisplay = displaySet(viewingDisplay,'wave',p.Results.wls);
viewingDisplay = displaySet(viewingDisplay,'spd primaries',displayGet(viewingDisplay,'spd primaries')*viewingDisplayScaleFactor);

% Sanity check by redoing the display correction based on the updates,
% utilizes the RenderAcrossDisplay functionality and ideally should be the
% same as before
imageRGB = gammaCorrection(imagergbLinear,startDisplay);
[imageRGBAcrossDisplays, ~ , ~] = ...
    RGBRenderAcrossDisplays(imageRGB, startDisplay, viewingDisplay, ...
    'viewingDisplayScaleFactor',p.Results.viewingDisplayScaleFactor, ...
    'linearInput',p.Results.linearInput,'verbose',p.Results.verbose, ...
    'scaleToMax',p.Results.scaleToMax,'SRGB',p.Results.SRGB, ...
    'wls', p.Results.wls);

% Select for a region that is half as large as the stimulus presentation
% region to be used when calculating EW, centered at the same location. The
% purpose of this is to account for any potential edge artifacts especially
% in recons by selecting for a central patch to calculate EW.
if isodd(length(idxXRange))
    centerPoint = ceil(length(idxXRange)/2);
    centerSpread = floor(length(idxXRange) / 4);
    idxXRangeCenter = idxXRange(centerPoint - centerSpread: centerPoint + centerSpread);
else
    centerPoint = ceil(length(idxXRange)/2);
    centerSpread = floor(length(idxXRange) / 4);stimEWFormer
    idxXRangeCenter = idxXRange(centerPoint - centerSpread + 1: centerPoint + centerSpread);
end

% Exclude the background from EW calculations by selecting for the region
% that overlaps with the full stimulus position. Then do the same for the
% center region.
RGBtoEWFull = imageRGB(idxXRange, idxXRange, :);
RGBtoEWCenter = imageRGB(idxXRangeCenter, idxXRangeCenter, :);

% Use the image regions established above to calculate the EW. Again, using
% the RGB values that correspond to those input to the mono display, NOT
% using the values that were returned after correcting to render across
% displays (that's just for visualizing).
[imageInfo.imageEWFull, ~] = ...
    RGBToEquivWavelength(RGBtoEWFull, startDisplay);
[imageInfo.imageEWCenter, ~] = ...
    RGBToEquivWavelength(RGBtoEWCenter, startDisplay);

imageInfo.meanEWFull = int64(mean(imageInfo.imageEWFull, 'all'));
imageInfo.meanEWCenter = int64(mean(imageInfo.imageEWCenter, 'all'));

% Visualize the newly rendered image if desired
if p.Results.showFigs
    figure()
    imshow(imageRGBAcrossDisplays)

    % Set the title to be the EW value determined (center or full)
    if p.Results.useFullEW
        title(sprintf('Equivalent Wavelength: %d', imageInfo.meanEWFull));
    else
        title(sprintf('Equivalent Wavelength: %d', imageInfo.meanEWCenter));
    end
end

% Store a version of the RGB image that was rendered across displays
imageInfo.imageRGBAcrossDisplays = imageRGBAcrossDisplays;

end