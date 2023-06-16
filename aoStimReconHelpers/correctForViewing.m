function cfv = correctForViewing(inputImageRGB, startDisplayName, ...
    viewingDisplayName, viewingDisplayScaleFactor, aoReconDir, ...
    displayGammaBits, displayGammaGamma, fieldSizeDegs, ...
    inputImageScaleFactor, idxXRange, idxYRange)
% Synopsis:
%    Correct input image to a form that better approximates how it would be
%    displayed on the start Display Monitor
%    inputImageRGB: RGB image matrix prior to gamma correction 
%    startDisplayName: 'mono' or 'conventional'
%    viewingDisplayName: 'mono' or 'conventional'
%    viewingDisplayScaleFactor: Some convenience scaling, value = 3 works
%    aoReconDir, displayGammaBits, displayGammaGamma: Params from pr struct
%
% Description:
%    Input a linear image specifying which monitor the image is being viewed under
%    (ViewingDisp) and apply adjustments based on which monitor the 
%    calculations were run under (startDisp). Ex: [1 1 0] RGB image appears
%    yellow under mono monitor for calculations but orangish under conventional. 
%
% See also: t_renderMonoDisplayImage, aoStimRecon, aoStimReconRunMany

% History:
%   04/04/23  chr  Made callable function from t_renderMonoDisplayImage

% 
% startDisplayName = 'mono';
% viewingDisplayName = 'conventional';

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

% Fix up gamma and wavelengths sampling if needed for starting display 
if (startOverwriteDisplayGamma)
    startGammaInput = linspace(0,1,2^displayGammaBits)';
    startGammaOutput = startGammaInput.^displayGammaGamma;
    startDisplay.gamma = startGammaOutput(:,[1 1 1]);
end

% Fix up gamma and wavelengths sampling if needed for viewing display
if (viewingOverwriteDisplayGamma)
    viewingGammaInput = linspace(0,1,2^displayGammaBits)';
    viewingGammaOutput = viewingGammaInput.^displayGammaGamma;
    viewingDisplay.gamma = viewingGammaOutput(:,[1 1 1]);
end

% Convenience assumption that all displays operate on the same wave structure
wls = (400:10:700)';
startDisplay = displaySet(startDisplay,'wave',wls);
viewingDisplay = displaySet(viewingDisplay,'wave',wls);

% Scale recon display primaries to try to keep things in range
viewingDisplay = displaySet(viewingDisplay,'spd primaries',displayGet(viewingDisplay,'spd primaries')*viewingDisplayScaleFactor);


meanLuminanceCdPerM2 = [];
[startScene, ~, startImageLinear] = sceneFromFile(inputImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, startDisplay);

%%%%%%%%
startImageRGB = gammaCorrection(startImageLinear * inputImageScaleFactor, startDisplay);
startScene = sceneSet(startScene, 'fov', fieldSizeDegs);
% [startScene, ~, startImageLinear] = sceneFromFile(startImageRGB, 'rgb', ...
%     meanLuminanceCdPerM2, startDisplay);
%%%%%%%%


% Get information we need to render scenes from their spectra through
% the recin display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mstart_rgbToXYZ = T_XYZ*displayGet(startDisplay,'spd primaries')*(wls(2)-wls(1));
Mstart_XYZTorgb = inv(Mstart_rgbToXYZ);
Mviewing_rgbToXYZ = T_XYZ*displayGet(viewingDisplay,'spd primaries')*(wls(2)-wls(1));
Mviewing_XYZTorgb = inv(Mviewing_rgbToXYZ);

% Render the scene from the forward display on the recon display, to try to
% match XYZ.
[thestartImagergbCalFormat,m,n] = ImageToCalFormat(startImageLinear);
thestartImageXYZCalFormat = Mstart_rgbToXYZ*thestartImagergbCalFormat;
theViewingImagergbCalFormat = Mviewing_XYZTorgb*thestartImageXYZCalFormat;
theViewingImagergb = CalFormatToImage(theViewingImagergbCalFormat,m,n);

% Truncate to keep within bounds 
minr = min(min(theViewingImagergb(:,:,1)));
ming = min(min(theViewingImagergb(:,:,2)));
minb = min(min(theViewingImagergb(:,:,3)));
maxr = max(max(theViewingImagergb(:,:,1)));
maxg = max(max(theViewingImagergb(:,:,2)));
maxb = max(max(theViewingImagergb(:,:,3)));
theViewingImagergbTruncated = theViewingImagergb;
theViewingImagergbTruncated(theViewingImagergbTruncated < 0) = 0;
theViewingImagergbTruncated(theViewingImagergbTruncated > 1) = 1;
sumBounds = [minr ming minb; maxr maxg maxb];

outputImageRGB = gammaCorrection(theViewingImagergbTruncated, viewingDisplay);


% boostScale = 1/max(theViewingImagergbTruncated, [], 'all');
% Boost the values based on the stimulus region (avoids being affected by
% random bright spots in the outer fringes)
boostScale = 1/max(theViewingImagergbTruncated(idxYRange, idxXRange, :), [], 'all');
theViewingImageBoosted = theViewingImagergbTruncated .* boostScale;

theViewingImageBoosted(theViewingImageBoosted < 0) = 0;
theViewingImageBoosted(theViewingImageBoosted > 1) = 1;
outputImageRGBBoost = gammaCorrection(theViewingImageBoosted,viewingDisplay); 


% Pull out the information for summary statistics using the INPUT
% uncorrected values
rgbStats = cell(3,4);
rgbStats(1,1) = {startImageLinear(idxYRange,idxXRange, 1)};
rgbStats(2,1) = {startImageLinear(idxYRange,idxXRange, 2)};
rgbStats(3,1) = {startImageLinear(idxYRange,idxXRange, 3)};

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


cfv.imageRGB             = outputImageRGB;
cfv.imageRGBBoost        = outputImageRGBBoost;
cfv.imageRGBBoostNoGamma = (theViewingImagergbTruncated .* boostScale);


cfv.imageRGBNoGamma      = theViewingImagergbTruncated;
cfv.bounds               = sumBounds;
cfv.rgbStats             = rgbStats;
cfv.viewingDisplay       = viewingDisplay; 


end