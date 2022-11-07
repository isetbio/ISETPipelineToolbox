function renderStructure = buildRenderStruct(aoReconDir, ...
    eccXDegs, eccYDegs, fieldSizeDegs, nPixels, pupilDiamMM, aoRender, ...
    defocusDiopters, overwriteDisplayGamma, displayName, displayFieldName, ...
    displayGammaBits, displayGammaGamma, randSeed, replaceCones, startCones, ...
    newCones, eccVars, subjectID, zernikeDataBase, quads)
% Synopsis:
%    Build render matrix if desired/needed
%
% Description:
%    Run this function if you would like to rebuild a new mosaic and 
%    render matrix.  This also gets run if there is no cached file corresponding
%    to the desired parameters. Once built, this file can be loaded from cache
%    for quicker running.
%
% See also: aoStimRecon, aoStimReconRunMany

% History:
%   08/16/22  chr  Made it a callable function
%   08/25/22  chr  Included portion for dichromacy
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options

% Get display
theDisplayLoad = load(fullfile(aoReconDir, 'displays', [displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' displayFieldName ';']);
if (overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^displayGammaBits)';
    gammaOutput = gammaInput.^displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end
clear theDisplayLoad;

% Spline underlying display wavelength down to the wavelength sampling we
% will eventually use in the calculations.
wls = (400:10:700)';
theDisplay = displaySet(theDisplay,'wave',wls);

% Create and setup cone mosaic
% 
% For AO, we put in subjectID == 0 which causes the zcoeffs to be all zero
% except for any specified defocus.
if (aoRender)
    if (eccVars)
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'wave', wls, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'subjectID', 0, ...
            'noLCA', true, ...
            'zernikeDataBase', zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', wls, ... 
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', 0, ...
            'noLCA', true, ...
            'zernikeDataBase', zernikeDataBase);
    end

% We build a normal optics structure. Allow specified defocus.
else
    if (eccVars)
        % Build normal optics structure. 
        theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'wave', wls, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, ...
            'subjectID', subjectID, ...
            'noLCA', false, ...
            'zernikeDataBase', zernikeDataBase);
    else
       theConeMosaic = ConeResponseCmosaic(eccXDegs, eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', wls, ... 
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', subjectID, ...
            'noLCA', false, ...
            'zernikeDataBase', zernikeDataBase);
    end
end

% Option to replace cones in mosaic with another kind to simulate
% dichromacy. 
if (replaceCones)
    if (quads(1).value)
        for q=2:5
            quads(q).regionCount = length(quads(q).percentL);
            quads(q).regionSpread = abs(0 - fieldSizeDegs / 2) / (quads(q).regionCount * 2 - 1);
            for r=1:quads(q).regionCount
                regionCones = find(...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > min(quads(q).xbounds) & ...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > min(quads(q).ybounds) & ...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < max(quads(q).xbounds) & ...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < max(quads(q).ybounds));
                
                indTracker = false(1,length(regionCones));
                newLAmount = round(length(regionCones) * quads(q).percentL(r));
                newLRegionInd = randperm(length(regionCones), newLAmount);
                indTracker(newLRegionInd) = true;
                
                newLMosaicInd = regionCones(indTracker);
                newMMosaicInd = regionCones(~indTracker);
                
                if ~isempty(newLMosaicInd)
                    theConeMosaic.Mosaic.reassignTypeOfCones(newLMosaicInd, cMosaic.LCONE_ID)
                end

                if ~isempty(newMMosaicInd)
                    theConeMosaic.Mosaic.reassignTypeOfCones(newMMosaicInd, cMosaic.MCONE_ID)
                end

                quads(q).xbounds = quads(q).xbounds + [quads(q).regionSpread -quads(q).regionSpread]; 
                quads(q).ybounds = quads(q).ybounds + [quads(q).regionSpread -quads(q).regionSpread];
            end
        end
    else
        for i=1:length(startCones)
            coneInd = find(theConeMosaic.Mosaic.coneTypes == startCones(i));
            theConeMosaic.Mosaic.reassignTypeOfCones(coneInd, newCones);
        end
    end
end

disp(['------------------------------------'])
% Generate render matrix
theConeMosaic.Display = theDisplay;
renderMatrix = theConeMosaic.forwardRender([nPixels nPixels 3], ...
    true, true, 'useDoublePrecision', true);
renderMatrix = double(renderMatrix);

% Push new info back into structure and save
renderStructure.theDisplay = theDisplay;
renderStructure.renderMatrix = renderMatrix;
renderStructure.theConeMosaic = theConeMosaic;
renderStructure.fieldSizeDegs = fieldSizeDegs;
renderStructure.eccX = eccXDegs;
renderStructure.eccY = eccYDegs;
renderStructure.nPixels = nPixels;
renderStructure.pupilDiamMM = pupilDiamMM;
renderStructure.AORender = aoRender;
renderStructure.defocusDiopters = defocusDiopters;
renderStructure.quadSeqInfo = quads;
end
    
