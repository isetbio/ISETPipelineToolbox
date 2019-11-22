function tonemapWrapper(imageSequence, figureNo)
if ~exist(figureNo, 'var')
    figureNo = 1;
end

tonemapParams.alpha = 0.1;
tonemapParams.gain = 70;

tonemapImageSet(figureNo, imageSequence, CRT12BitDisplay, tonemapParams);
end


function tonemapImageSet(figureNo, imageSetRGB, display, tonemapParams)


% Retrieve display's XYZ -> RGB matric
displayXYZtoRGBMatrix = inv(displayGet(display, 'rgb2xyz'));


% Allocate memory
imageNums = size(imageSetRGB,1);
mRows = size(imageSetRGB,2);
nCols = size(imageSetRGB,3);
imageSetLuma   = zeros(imageNums, mRows, nCols);
imageSetChroma = zeros(imageNums, mRows, nCols, 2);


% Decompose images into luma/chroma channels
setLuminances = [];
for imageNo = 1:imageNums
    [imageSetLuma(imageNo,:,:), imageSetChroma(imageNo,:,:,:)] = ...
        RGBtoLumaChroma(squeeze(imageSetRGB(imageNo,:,:,:)),displayXYZtoRGBMatrix);
    
    
    setLuminances = cat(1, setLuminances, imageSetLuma(imageNo,:,:));
end


% Compute all scenes key
setKey = computeSetKey(setLuminances(:));


for imageNo = 1:imageNums
    % Tone map the luminance channel
    imageSetToneMappedLuma(imageNo,:,:) = toneMapLuma(squeeze(imageSetLuma(imageNo,:,:)), setKey, tonemapParams);
    
    
    % Recostuct RGB image using tonemapped luminance channel
    [imageSetRGBToneMapped(imageNo,:,:,:), outOfGamutImage1, outOfGamutImage2] = primariesFromChromaticitiesAndLuminances(...
        squeeze(imageSetToneMappedLuma(imageNo,:,:)), ...
        squeeze(imageSetChroma(imageNo,:,:,:)), ...
        displayXYZtoRGBMatrix);
    
    
    % Display results
    displayResults(figureNo, imageNo, imageNums, ...
        squeeze(imageSetRGB(imageNo,:,:,:)), ...
        squeeze(imageSetRGBToneMapped(imageNo,:,:,:)), ...
        squeeze(imageSetLuma(imageNo,:,:)), ...
        squeeze(imageSetRGBToneMapped(imageNo,:,:,2)), ...
        outOfGamutImage1, outOfGamutImage2);
end
end


function displayResults(figureNo, imageNo, imageNums, imageRGB, imageRGBToneMapped, imageLuma, imageToneMappedLuma, outOfGamutImage1, outOfGamutImage2)


if (imageNo == 1)
    hFig = figure(figureNo); clf;
    set(hFig, 'Position', [10 10 890 1500]);
end


subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', 6, ...
    'colsNum', imageNums, ...
    'heightMargin',  0.04, ...
    'widthMargin',    0.03, ...
    'leftMargin',     0.04, ...
    'rightMargin',    0.01, ...
    'bottomMargin',   0.03, ...
    'topMargin',      0.01);


subplot('Position',subplotPosVectors(1,imageNo).v);
imshow(imageRGB, [0 1]);
title('original image')


subplot('Position',subplotPosVectors(2,imageNo).v);
imshow(imageRGBToneMapped, [0 1]);
title('tone mapped image')


subplot('Position',subplotPosVectors(3,imageNo).v);
plot(imageLuma(:), imageToneMappedLuma(:), 'k.')
axis 'square';
xlabel('input luminance');
ylabel('tone mapped luminance');
title('luminance mapping');


subplot('Position',subplotPosVectors(4,imageNo).v);
inR = squeeze(imageRGB(:,:,1));
outR = squeeze(imageRGBToneMapped(:,:,1));
plot(inR(:), outR(:), 'r.'); hold on;


inG= squeeze(imageRGB(:,:,2));
outG = squeeze(imageRGBToneMapped(:,:,2));
plot(inG(:), outG(:), 'g.');


inB= squeeze(imageRGB(:,:,3));
outB = squeeze(imageRGBToneMapped(:,:,3));
plot(inB(:), outB(:), 'b.');


axis 'square';
set(gca, 'XLim', [0 1], 'YLim', [0 1]);
xlabel('input primary');
ylabel('output primary');
title('RGB mapping');


subplot('Position',subplotPosVectors(5,imageNo).v);
imshow(outOfGamutImage1);
title(sprintf('out of gamut, clipped to 0 (%d pixels)', numel(find(outOfGamutImage1(:) > 0))));


subplot('Position',subplotPosVectors(6,imageNo).v);
imshow(outOfGamutImage2);
title(sprintf('out of gamut, clipped to 1 (%d pixels)', numel(find(outOfGamutImage2(:) > 0))));


end




function [RGBprimaries, OutOfGamutImage1, OutOfGamutImage2] = primariesFromChromaticitiesAndLuminances(luminances, chromaticities, displayXYZtoRGBMatrix)
[mRows, nCols] = size(luminances);
xyY = cat(2, reshape(chromaticities, [nCols*mRows 2]), reshape(luminances, [nCols*mRows 1]));
XYZCalFormat = xyYToXYZ(xyY');
RGBprimariesCalFormat = displayXYZtoRGBMatrix * XYZCalFormat;
RGBprimaries = CalFormatToImage(RGBprimariesCalFormat, nCols, mRows);
OutOfGamutImage1 = zeros(mRows, nCols);
OutOfGamutImage2 = zeros(mRows, nCols);
for channelIndex = 1:size(RGBprimaries,3)
    theChannelData = squeeze(RGBprimaries(:,:,channelIndex));
    idx = find(theChannelData(:)<0);
    OutOfGamutImage1(idx) = 1;
    idx = find(theChannelData(:)>1);
    OutOfGamutImage2(idx) = 1;
end
RGBprimaries(RGBprimaries < 0) = 0;
RGBprimaries(RGBprimaries > 1) = 1;
end


function toneMappedImageLuminance = toneMapLuma(imageLuminance, setKey, tonemapParams)
% Deal with negative luminances. Truncate to 0
imageLuminance(imageLuminance<0) = 0;


% Reinhard tone mapping
% Scale luminance according to alpha parameter and scene key
scaledInputLuminance = tonemapParams.alpha / setKey * imageLuminance;
% Compress high luminances
toneMappedImageLuminance = tonemapParams.gain * scaledInputLuminance ./ (1.0+scaledInputLuminance);
end


function key = computeSetKey(setLuminances)
setLuminances(setLuminances<0) = 0;
delta = 0.0001; % small delta to avoid taking log(0) when encountering pixels with zero luminance
key = exp((1/numel(setLuminances))*sum(log(setLuminances + delta)));
end


function [imageLuma, imageChroma] = RGBtoLumaChroma(rgbPrimariesImage,displayXYZtoRGBMatrix)


% Transform the RGB image [Rows x Cols x 3] into a [3 x N] matrix for faster computations
[rgbPrimariesCalFormat,nCols,mRows] = ImageToCalFormat(rgbPrimariesImage);


XYZCalFormat = inv(displayXYZtoRGBMatrix) * rgbPrimariesCalFormat;


% xyY values of the image
xyYCalFormat = XYZToxyY(XYZCalFormat);
imageChroma = CalFormatToImage(xyYCalFormat(1:2,:), nCols,mRows);
imageLuma = CalFormatToImage(xyYCalFormat(3,:),  nCols,mRows);
end