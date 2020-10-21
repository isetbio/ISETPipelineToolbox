function rgbImage = invGammaCorrection(linearImage, display)
%INVGAMMACORRECTION Summary of this function goes here
%   Detailed explanation goes here

gammaTable = displayGet(display, 'gamma table');
nInputLevels = size(gammaTable,1);
gammaInput   = linspace(0,1,nInputLevels);

PTBcal = ptb.GeneratePsychToolboxCalStruct(...
                'name', displayGet(display, 'name'), ...
                'gammaInput', gammaInput, ...
                'gammaTable', gammaTable, ...
                'wave', displayGet(display, 'wave'), ...
                'spd',  displayGet(display, 'spd'), ...
                'ambientSpd', zeros(length(displayGet(display, 'wave')),1));

gammaMethod = 1;
PTBcal = SetGammaMethod(PTBcal, gammaMethod, nInputLevels);
[linearizedCalFormat,m,n] = ImageToCalFormat(linearImage);
rgbImage = CalFormatToImage(PrimaryToSettings(PTBcal,linearizedCalFormat),m,n);

end