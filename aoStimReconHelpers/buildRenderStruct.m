function renderStructure = buildRenderStruct(pr, cnv, pupilDiamMM, ...
    aoRender, noLCA, defocusDiopters, randSeed, replaceCones, startCones, ...
    newCones, eccVars, subjectID, zernikeDataBase, chrom)
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
theDisplayLoad = load(fullfile(pr.aoReconDir, 'displays', [pr.displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' cnv.displayFieldName ';']);
if (cnv.overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^pr.displayGammaBits)';
    gammaOutput = gammaInput.^pr.displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end
clear theDisplayLoad;

% Spline underlying display wavelength down to the wavelength sampling we
% will eventually use in the calculations.
wls = (400:10:700)';
theDisplay = displaySet(theDisplay,'wave',wls);


fieldSizeDegs = pr.fieldSizeMinutes/60;

% Create and setup cone mosaic
% 
% For AO, we put in subjectID == 0 which causes the zcoeffs to be all zero
% except for any specified defocus.
if (aoRender)
    if (eccVars)
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', wls, ... 
            'subjectID', 0, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
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
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    end

% We build a normal optics structure. Allow specified defocus.
else
    if (eccVars)
        % Build normal optics structure. 
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', pupilDiamMM, 'useRandomSeed', randSeed, ...
            'defocusDiopters',defocusDiopters, 'wave', wls, ... 
            'subjectID', subjectID, ...
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    else
       theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
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
            'noLCA', noLCA, ...
            'zernikeDataBase', zernikeDataBase);
    end
end

% Option to replace cones in mosaic with another kind to simulate
% dichromacy. 
if (replaceCones)
    if (pr.quads(6).value)
        storedDir = fullfile(pr.aoReconDir, pr.versEditor);
        theConeMosaic = overrideQuads(theConeMosaic, pr.eccXDegs, pr.eccYDegs, chrom, storedDir);
    elseif (pr.quads(1).value)
        for q=2:5
            % Ensure input values are equivalent
            if length(pr.quads(q).percentL) ~= length(pr.quads(q).percentS)
                error(['Unequal number of regions across cone classes in ', ...
                    pr.quads(q).name])
            end

            % Collect the number of regions and equidistant sizes from
            % number of L percentage values entered
            pr.quads(q).regionCount = length(pr.quads(q).percentL);
            pr.quads(q).regionSpread = abs(0 - fieldSizeDegs / 2) / (pr.quads(q).regionCount * 2 - 1);
            for r=1:pr.quads(q).regionCount
                % Capture the indices of cones that fall within each region
                regionCones = find(...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > min(pr.quads(q).xbounds) & ...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > min(pr.quads(q).ybounds) & ...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < max(pr.quads(q).xbounds) & ...
                theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < max(pr.quads(q).ybounds));
                
                % Initialize a index tracker using boolean vectors for
                % efficiency in L/M random assignments
                indTracker = false(1,length(regionCones));

                % Convert the desired L percentage to a number of cones
                newLAmount = round(length(regionCones) * pr.quads(q).percentL(r));

                % Select the desired number of cones from the region to be
                % converted to L cones, represented as true in the tracker
                newLRegionInd = randperm(length(regionCones), newLAmount);
                indTracker(newLRegionInd) = true;

                % Convert the desired S percentage to a number of cones and
                % apply this number to the mosaic regions
                newSAmount = round(length(regionCones) * pr.quads(q).percentS(r));
                newSRegionInd = randperm(length(regionCones), newSAmount);

                % Apply desired percentages to the L and S cones using the
                % tracker. Then override all with aplied S cone values
                newLMosaicInd = regionCones(indTracker);
                newMMosaicInd = regionCones(~indTracker);
                newSMosaicInd = regionCones(newSRegionInd);

                % If cone types should be present, complete the switch
                if ~isempty(newLMosaicInd)
                    theConeMosaic.Mosaic.reassignTypeOfCones(newLMosaicInd, cMosaic.LCONE_ID)
                end
                if ~isempty(newMMosaicInd)
                    theConeMosaic.Mosaic.reassignTypeOfCones(newMMosaicInd, cMosaic.MCONE_ID)
                end
                if ~isempty(newSMosaicInd)
                    theConeMosaic.Mosaic.reassignTypeOfCones(newSMosaicInd, cMosaic.SCONE_ID)
                end
                
                % Update the region boundaries to move inward
                pr.quads(q).xbounds = pr.quads(q).xbounds + [pr.quads(q).regionSpread -pr.quads(q).regionSpread]; 
                pr.quads(q).ybounds = pr.quads(q).ybounds + [pr.quads(q).regionSpread -pr.quads(q).regionSpread];
            end
        end
    else
        % Preset mosaics for more typical conversions. 
        for i=1:length(startCones)
            coneInd = find(theConeMosaic.Mosaic.coneTypes == startCones(i));
            theConeMosaic.Mosaic.reassignTypeOfCones(coneInd, newCones);
        end
    end
end

% Generate render matrix
theConeMosaic.Display = theDisplay;
renderMatrix = theConeMosaic.forwardRender([pr.nPixels pr.nPixels 3], ...
    true, true, 'useDoublePrecision', true);
renderMatrix = double(renderMatrix);

% Push new info back into structure and save
renderStructure.theDisplay = theDisplay;
renderStructure.renderMatrix = renderMatrix;
renderStructure.theConeMosaic = theConeMosaic;
renderStructure.fieldSizeDegs = fieldSizeDegs;
renderStructure.eccX = pr.eccXDegs;
renderStructure.eccY = pr.eccYDegs;
renderStructure.nPixels = pr.nPixels;
renderStructure.pupilDiamMM = pupilDiamMM;
renderStructure.AORender = aoRender;
renderStructure.defocusDiopters = defocusDiopters;
renderStructure.quadSeqInfo = pr.quads;
end
    
