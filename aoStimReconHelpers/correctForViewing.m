function cfv = correctForViewing(inputImagergb, startDisplayName, ...
    viewingDisplayName, viewingDisplayScaleFactor, aoReconDir, ...
    displayGammaBits, displayGammaGamma, fieldSizeDegs, ...
    inputImageScaleFactor, idxXRange, idxYRange)
% Correct input image for better display
%
% Synopsis:
%     cfv = correctForViewing(inputImagergb, startDisplayName, ...
%         viewingDisplayName, viewingDisplayScaleFactor, aoReconDir, ...
%         displayGammaBits, displayGammaGamma, fieldSizeDegs, ...
%         inputImageScaleFactor, idxXRange, idxYRange)
%
% Description:
%    Input a linear image specifying which monitor the image is being viewed under
%    (ViewingDisp) and apply adjustments based on which monitor the 
%    calculations were run under (startDisp). Ex: [1 1 0] RGB image appears
%    yellow under mono monitor for calculations but orangish under conventional.
%
% Inputs:
%    inputImagergb:             RGB image matrix prior to gamma correction 
%    startDisplayName:          'mono' or 'conventional'
%    viewingDisplayName:        'mono' or 'conventional'
%    viewingDisplayScaleFactor: Some convenience scaling, value = 3 works
%    aoReconDir:
%    displayGammaBits:
%    displayGammaGamma:
%    fieldSizeDegs:             Does not appear to be used
%    inputImageScaleFactor:     Does not appear to be used
%    idxXRange:                 Used to define region under stimulus in input image
%    idxYRange:                 Used to define region under stimulus in input image
%
% Outputs:
%   cfv:                        A structure that Carlos understands                   
%
% See also: t_renderMonoDisplayImage, RGBRenderAcrossDisplays, aoStimRecon, aoStimReconRunMany

% History:
%   04/04/23  chr  Made callable function from t_renderMonoDisplayImage
%   10/12/23  dhb  Call through RGBRenderAcrossDisplays
%   03/29/24  chr  THIS FUNCTION IS NO LONGER ACTIVE, STILL HERE TO
%                  ENSURE LOSING IT DOESN'T BREAK THINGS, BUT PLAN IS TO 
%                  DELETE AT EARLIEST CONVENIENCE

% Initialize the struct for holding cfv
cfv = struct; 

% Specify variables depending on the starting display 
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

% Get displays
startDisplayLoad = load(fullfile(aoReconDir, 'displays', [startDisplayName 'Display.mat']));
eval(['startDisplay = startDisplayLoad.' startDisplayFieldName ';']);
viewingDisplayLoad = load(fullfile(aoReconDir, 'displays', [viewingDisplayName 'Display.mat']));
eval(['viewingDisplay = viewingDisplayLoad.' viewingDisplayFieldName ';']);
clear startDisplayLoad viewingDisplayLoad

% Fix up gamma and wavelength sampling if needed for starting display 
if (startOverwriteDisplayGamma)
    startGammaInput = linspace(0,1,2^displayGammaBits)';
    startGammaOutput = startGammaInput.^displayGammaGamma;
    startDisplay.gamma = startGammaOutput(:,[1 1 1]);
end

% Fix up gamma and wavelength sampling if needed for viewing display
if (viewingOverwriteDisplayGamma)
    viewingGammaInput = linspace(0,1,2^displayGammaBits)';
    viewingGammaOutput = viewingGammaInput.^displayGammaGamma;
    viewingDisplay.gamma = viewingGammaOutput(:,[1 1 1]);
end

% Call underlying function to do the work. Can use SRGB or not
% !!! Note, as it is written now, will not work because passing in a linear
% rgb value when fxn requires gamma corrected RGB value. 
SRGB = true;
[outputImageRGB,theViewingImagergbTruncated,theViewingImagergb] = RGBRenderAcrossDisplays(inputImagergb, startDisplay, [], ...
            'viewingDisplayScaleFactor',viewingDisplayScaleFactor, ...
            'linearInput',true,'verbose',false, ...
            'scaleToMax',false,'SRGB',SRGB);

% Collect stats on the untruncted version
minr = min(min(theViewingImagergb(:,:,1)));
ming = min(min(theViewingImagergb(:,:,2)));
minb = min(min(theViewingImagergb(:,:,3)));
maxr = max(max(theViewingImagergb(:,:,1)));
maxg = max(max(theViewingImagergb(:,:,2)));
maxb = max(max(theViewingImagergb(:,:,3)));
sumBounds = [minr ming minb; maxr maxg maxb];

% Scale the values based on the stimulus region (avoids being affected by
% random bright spots in the outer fringes)
boostScale = 1/max(theViewingImagergbTruncated(idxYRange, idxXRange, :), [], 'all');
theViewingImagergbBoosted = theViewingImagergbTruncated .* boostScale;
theViewingImagergbBoosted(theViewingImagergbBoosted > 1) = 1;
if (SRGB)
    outputImageRGBBoost  = double(SRGBGammaCorrect(theViewingImagergbBoosted,0))/255;
else
    outputImageRGBBoost = gammaCorrection(theViewingImagergbBoosted, viewingDisplay);
end

% Pull out the information for summary statistics using the input
% uncorrected linear rgb values
rgbStats = cell(3,5);
rgbStats(1,1) = {inputImagergb(idxYRange,idxXRange, 1)};
rgbStats(2,1) = {inputImagergb(idxYRange,idxXRange, 2)};
rgbStats(3,1) = {inputImagergb(idxYRange,idxXRange, 3)};

rgbStats(1,2) = {mean(rgbStats{1,1}(:))};
rgbStats(2,2) = {mean(rgbStats{2,1}(:))};
rgbStats(3,2) = {mean(rgbStats{3,1}(:))};

rgbStats(1,3) = {std(rgbStats{1,1}(:))};
rgbStats(2,3) = {std(rgbStats{2,1}(:))};
rgbStats(3,3) = {std(rgbStats{3,1}(:))};

rgbStats(1,4) = {ones(size(rgbStats{1,1})) * rgbStats{1,2}};
rgbStats(2,4) = {ones(size(rgbStats{2,1})) * rgbStats{2,2}};
rgbStats(3,4) = {ones(size(rgbStats{3,1})) * rgbStats{3,2}};

rgbStats(1,5) = {rgbStats{1,2} / (rgbStats{1,2} + rgbStats{2,2})};

% Set up return structure
cfv.imageRGB             = outputImageRGB;
cfv.imageRGBBoost        = outputImageRGBBoost;
cfv.imageRGBBoostNoGamma = (theViewingImagergbTruncated .* boostScale);
cfv.imageRGBNoGamma      = theViewingImagergbTruncated;
cfv.bounds               = sumBounds;
cfv.rgbStats             = rgbStats;
cfv.viewingDisplay       = viewingDisplay; 

end