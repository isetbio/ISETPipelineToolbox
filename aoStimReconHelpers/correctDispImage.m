function [outputImageRGB, sumBounds, outputImageRGBBoost] = ...
    correctDispImage(inputImageRGB, trueDisplayName, ...
    viewingDisplayName, viewingDisplayScaleFactor, aoReconDir, ...
    displayGammaBits, displayGammaGamma, fieldSizeDegs, inputImageScaleFactor)
% Synopsis:
%    Correct input image to a form that better approximates how it would be
%    displayed on the true Display Monitor
%    inputImageRGB: RGB image matrix prior to gamma correction 
%    trueDisplayName: 'mono' or 'conventional'
%    viewingDisplayName: 'mono' or 'conventional'
%    viewingDisplayScaleFactor: Some convenience scaling, value = 3 works
%    aoReconDir, displayGammaBits, displayGammaGamma: Params from pr struct
%
% Description:
%    Input a linear image specifying which monitor the image is being viewed under
%    (ViewingDisp) and apply adjustments based on which monitor the 
%    calculations were run under (TrueDisp). Ex: [1 1 0] RGB image appears
%    yellow under mono monitor for calculations but orangish under conventional. 
%
% See also: t_renderMonoDisplayImage, aoStimRecon, aoStimReconRunMany

% History:
%   04/04/23  chr  Made callable function from t_renderMonoDisplayImage

% 
% trueDisplayName = 'mono';
% viewingDisplayName = 'conventional';


% Set up display specific fields
switch (trueDisplayName)
    case 'conventional'
        trueDisplayFieldName = 'CRT12BitDisplay';
        trueOverwriteDisplayGamma = false;
    case 'mono'
        trueDisplayFieldName = 'monoDisplay';
        trueOverwriteDisplayGamma = true;
    otherwise
        error('Unknown recon display specified');
end
switch (viewingDisplayName)
    case 'conventional'
        viewingDisplayFieldName = 'CRT12BitDisplay';
        viewingOverwriteDisplayGamma = false;
    case 'mono'
        viewingDispFieldName = 'monoDisplay';
        viewingOverwriteDisplayGamma = true;
    otherwise
        error('Unknown forward display specified');
end


% Get displays
trueDisplayLoad = load(fullfile(aoReconDir, 'displays', [trueDisplayName 'Display.mat']));
eval(['trueDisplay = trueDisplayLoad.' trueDisplayFieldName ';']);
viewingDisplayLoad = load(fullfile(aoReconDir, 'displays', [viewingDisplayName 'Display.mat']));
eval(['viewingDisplay = viewingDisplayLoad.' viewingDisplayFieldName ';']);
clear trueDisplayLoad viewingDisplayLoad

% Fix up gamma and wavelengths sampling
if (trueOverwriteDisplayGamma)
    trueGammaInput = linspace(0,1,2^displayGammaBits)';
    trueGammaOutput = trueGammaInput.^displayGammaGamma;
    trueDisplay.gamma = trueGammaOutput(:,[1 1 1]);
end

if (viewingOverwriteDisplayGamma)
    viewingGammaInput = linspace(0,1,2^displayGammaBits)';
    viewingGammaOutput = viewingGammaInput.^displayGammaGamma;
    viewingDisplay.gamma = viewingGammaOutput(:,[1 1 1]);
end

% Convenience assumption that all displays operate on the same wave structure
wls = (400:10:700)';
trueDisplay = displaySet(trueDisplay,'wave',wls);
viewingDisplay = displaySet(viewingDisplay,'wave',wls);

% Scale recon display primaries to try to keep things in range
viewingDisplay = displaySet(viewingDisplay,'spd primaries',displayGet(viewingDisplay,'spd primaries')*viewingDisplayScaleFactor);


meanLuminanceCdPerM2 = [];
[trueScene, ~, trueImageLinear] = sceneFromFile(inputImageRGB, 'rgb', ...
    meanLuminanceCdPerM2, trueDisplay);
trueImageRGB = gammaCorrection(trueImageLinear * inputImageScaleFactor, trueDisplay);
trueScene = sceneSet(trueScene, 'fov', fieldSizeDegs);

% [trueScene, ~, trueImageLinear] = sceneFromFile(trueImageRGB, 'rgb', ...
%     meanLuminanceCdPerM2, trueDisplay);


% Get information we need to render scenes from their spectra through
% the recin display.
theXYZStruct = load('T_xyz1931');
T_XYZ = SplineCmf(theXYZStruct.S_xyz1931,683*theXYZStruct.T_xyz1931,wls);
Mtrue_rgbToXYZ = T_XYZ*displayGet(trueDisplay,'spd primaries')*(wls(2)-wls(1));
Mtrue_XYZTorgb = inv(Mtrue_rgbToXYZ);
Mviewing_rgbToXYZ = T_XYZ*displayGet(viewingDisplay,'spd primaries')*(wls(2)-wls(1));
Mviewing_XYZTorgb = inv(Mviewing_rgbToXYZ);

% Render the scene from the forward display on the recon display, to try to
% match XYZ.
[theTrueImagergbCalFormat,m,n] = ImageToCalFormat(trueImageLinear);
theTrueImageXYZCalFormat = Mtrue_rgbToXYZ*theTrueImagergbCalFormat;
theViewingImagergbCalFormat = Mviewing_XYZTorgb*theTrueImageXYZCalFormat;
theViewingImagergb = CalFormatToImage(theViewingImagergbCalFormat,m,n);

% Truncate
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

boostScale = 1/max(theViewingImagergbTruncated, [], 'all');
outputImageRGBBoost = gammaCorrection(...
    (theViewingImagergbTruncated .* boostScale),viewingDisplay); 

% figure; clf; imshow(outputImageRGB);
% title(['Viewing Image'])
% keyboard

end