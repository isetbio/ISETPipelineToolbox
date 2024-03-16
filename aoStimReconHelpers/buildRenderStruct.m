function renderStructure = buildRenderStruct(pr, cnv, stage, varargin)
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
%% Set the stage
%
% Unpack the variables that change dependent on whether building a forward
% mosaic or recon mosaic
st = unpackStage(pr, cnv, stage);

%% Build base mosaic using input parameters
%
% Get display
theDisplayLoad = load(fullfile(pr.aoReconDir, 'displays', [pr.displayName 'Display.mat']));
eval(['theDisplay = theDisplayLoad.' cnv.displayFieldName ';']);
if (cnv.overwriteDisplayGamma)
    gammaInput = linspace(0,1,2^pr.displayGammaBits)';
    gammaOutput = gammaInput.^pr.displayGammaGamma;
    theDisplay.gamma = gammaOutput(:,[1 1 1]);
end
clear theDisplayLoad;

% Set some pertinent variables
wls = (400:1:700)';
theDisplay = displaySet(theDisplay,'wave',wls);
fieldSizeDegs = pr.fieldSizeMinutes / 60;

% Create and setup cone mosaic
%
% For AO, we put in subjectID == 0 which causes the zcoeffs to be all zero
% except for any specified defocus.
if (st.aoRender)
    if (st.eccVars)
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', st.pupilDiamMM, 'useRandomSeed', st.randSeed, ...
            'defocusDiopters',st.defocusDiopters, 'wave', wls, ...
            'subjectID', 0, ...
            'noLCA', st.noLCA, ...
            'zernikeDataBase', st.zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', st.pupilDiamMM, 'useRandomSeed', st.randSeed, ...
            'defocusDiopters',st.defocusDiopters, 'wave', wls, ...
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', 0, ...
            'noLCA', st.noLCA, ...
            'zernikeDataBase', st.zernikeDataBase);
    end

    % We build a normal optics structure. Allow specified defocus.
else
    if (st.eccVars)
        % Build normal optics structure.
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', st.pupilDiamMM, 'useRandomSeed', st.randSeed, ...
            'defocusDiopters',st.defocusDiopters, 'wave', wls, ...
            'subjectID', st.subjectID, ...
            'noLCA', st.noLCA, ...
            'zernikeDataBase', st.zernikeDataBase);
    else
        theConeMosaic = ConeResponseCmosaic(pr.eccXDegs, pr.eccYDegs, ...
            'fovealDegree', fieldSizeDegs, 'pupilSize', st.pupilDiamMM, 'useRandomSeed', st.randSeed, ...
            'defocusDiopters',st.defocusDiopters, 'wave', wls, ...
            'rodIntrusionAdjustedConeAperture', false, ...
            'eccVaryingConeAperture', false, ...
            'eccVaryingConeBlur', false, ...
            'eccVaryingOuterSegmentLength', false, ...
            'eccVaryingMacularPigmentDensity', false, ...
            'eccVaryingMacularPigmentDensityDynamic', false, ...
            'anchorAllEccVaryingParamsToTheirFovealValues', true, ...
            'subjectID', st.subjectID, ...
            'noLCA', st.noLCA, ...
            'zernikeDataBase', st.zernikeDataBase);
    end
end

%% Build render matrix from base mosaic
if ~pr.useCustomMosaic
    % Set filler mosaicConeInfo variable for now, should expand it to
    % actually return values as if mosaic was global
    mosaicConeInfo = [];

    % Build the render structure for the base mosaic
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
    renderStructure.pupilDiamMM = st.pupilDiamMM;
    renderStructure.AORender = st.aoRender;
    renderStructure.defocusDiopters = st.defocusDiopters;
    renderStructure.mosaicConeInfo = mosaicConeInfo;

    % Save it to the output location 
    save(fullfile(cnv.renderDir, 'xRenderStructures', st.renderStructureName),'renderStructure','-v7.3');
else
    %% Build render matrix based on edited mosaic
    % 
    % Establish dimensions of the full montage
    allRows = length(pr.stimSizeDegsDomain) * length(pr.focalRegionDomain);
    allColms = length(pr.focalPropLListDomain);

    for h = 1:length(pr.focalVariantDomain)
        for i = 1:length(pr.stimSizeDegsDomain)
            for j = 1:length(pr.focalRegionDomain)
                for k = 1:length(pr.focalPropLListDomain)
                    % Build the custom mosaics using setConeProportions
                    [theConeMosaic, mosaicConeInfo] = setConeProportions(pr.focalRegionDomain(j), ...
                        pr.focalPropLListDomain(k), pr.focalVariantDomain(h), theConeMosaic, pr.eccXDegs, pr.eccYDegs, ...
                        pr.stimSizeDegsDomain(i), pr.fieldSizeMinutes, pr.regionVariant, pr.propL, pr.propS);

                    % Build th render structure for the custom mosaic 
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
                    renderStructure.pupilDiamMM = st.pupilDiamMM;
                    renderStructure.AORender = st.aoRender;
                    renderStructure.defocusDiopters = st.defocusDiopters;
                    renderStructure.mosaicConeInfo = mosaicConeInfo;

                    % Save it to the output location
                    save(fullfile(cnv.renderDir, 'xRenderStructures', st.renderStructureName),'renderStructure','-v7.3');
                end
            end
        end
    end
end
end



