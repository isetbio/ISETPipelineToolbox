function renderStructure = buildRenderStruct(pr, cnv, pupilDiamMM, ...
    aoRender, noLCA, defocusDiopters, randSeed, replaceCones, startCones, ...
    newCones, eccVars, subjectID, zernikeDataBase, chrom, varargin)
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

p = inputParser;
p.addParameter('viewBounds', false, @islogical);
parse(p, varargin{:});

%% Variable check

% prVars = (["one" "two" "three" "four"]);
% prVarsCheck = isfield(pr, prVars);
% prVarsMissing = prVars(~prVarsCheck);
% if ~isempty(prVarsMissing)
%     error('Missing variables from pr struct')
% end
%
% cnvVars = (["one" "two" "three" "four"]);
% cnvVarsCheck = isfield(cnv, cnvVars);
% cnvVarsMissing = cnvVars(~cnvVarsCheck);
% if ~isempty(cnvVarsMissing)
%     error('Missing variables from cnv struct')
% end

%% Build Mosaic
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

% Spline underlying display wavelength down to the wavelength sampling we
% will eventually use in the calculations.
% This is at 10 but we're transitioning to 1 so maybe adjust here as well
% and see how that goes?
wls = (400:1:700)';
theDisplay = displaySet(theDisplay,'wave',wls);

fieldSizeDegs = pr.fieldSizeMinutes / 60;

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


%% Organize Variables
for w = 1:length(pr.focalRegionDomain)
    switch pr.focalRegionDomain(w)
        case "center"
            focalRegionPlot(w) = "Center";
        case "nearSurround"
            focalRegionPlot(w) = "Near";
        case "distantSurround"
            focalRegionPlot(w) = "Distant";
        otherwise
            focalRegionPlot(w) = pr.focalRegionDomain(w);
    end
end

allRows = length(pr.stimSizeDegsDomain) * length(pr.focalRegionDomain);
allColms = length(pr.focalPropLListDomain);

%% Build render matrix or view mosaic monpr.useBaseMosaic
    for h = 1:length(pr.focalVariantDomain)
        for i = 1:length(pr.stimSizeDegsDomain)
            for j = 1:length(pr.focalRegionDomain)
                for k = 1:length(pr.focalPropLListDomain)

                    [theConeMosaic, mosaicConeInfo] = setConeProportions(pr.focalRegionDomain(j), ...
                        pr.focalPropLListDomain(k), pr.focalVariantDomain(h), theConeMosaic, pr.eccXDegs, pr.eccYDegs, ...
                        pr.stimSizeDegsDomain(i), pr.fieldSizeMinutes);

                    % Build the Render Struct for the created Mosaic
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
                    renderStructure.mosaicConeInfo = mosaicConeInfo;
                    % renderStructure.quadSeqInfo = pr.quads;

                    save(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'renderStructure','-v7.3');
                end


            end
        end
    end
end
end


