function buildMosaicMontage(pr, cnv, stage, varargin)
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
%   03/13/24  chr  Split from buildRenderStruct

%% Set the stage
%
% Unpack the variables that change dependent on whether building a forward
% mosaic or recon mosaic
st = unpackStage(pr, cnv, stage);

%% DAVID THINKS THIS DISPLAY STUFF IS NOT NEEDED HERE
%% Build base mosaic using input parameters
%
% Get display
% theDisplayLoad = load(fullfile(pr.aoReconDir, 'displays', [pr.displayName 'Display.mat']));
% eval(['theDisplay = theDisplayLoad.' cnv.displayFieldName ';']);
% if (cnv.overwriteDisplayGamma)
%     gammaInput = linspace(0,1,2^pr.displayGammaBits)';
%     gammaOutput = gammaInput.^pr.displayGammaGamma;
%     theDisplay.gamma = gammaOutput(:,[1 1 1]);
% end
% clear theDisplayLoad;
% 
% % Set some pertinent variables
% wls = (400:1:700)';
% theDisplay = displaySet(theDisplay,'wave',wls);
fieldSizeDegs = pr.fieldSizeMinutes / 60;

% Create and setup base cone mosaic
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

%% Axes Naming
% 
% Adjust how variables are presented to make montage prettier
for w = 1:length(pr.focalRegionDomain)
    switch pr.focalRegionDomain(w)
        case "center"
            regionAxesNames(w) = "Center";
        case "nearSurround"
            regionAxesNames(w) = "Near";
        case "distantSurround"
            regionAxesNames(w) = "Distant";
        case "multiple"
            regionAxesNames(w) = "Multiple";
        case "global"
            regionAxesNames(w) = "Global";
    end
end

switch stage
    case "forward"
        stageTitle = "Forward";
    case "recon"
        stageTitle = "Recon";
end

%% Build mosaic montage based on edited mosaic
%
% Establish dimensions of the full montage
allRows = length(pr.stimSizeDegsDomain) * length(pr.focalRegionDomain);
allColms = length(pr.focalPropLListDomain);
viewBounds = false; %%%%%%%%%%%%%
if pr.useCustomMosaic
    for h = 1:length(pr.focalVariantDomain)

        theFig = figure;
        t = tiledlayout(allRows, allColms, 'TileSpacing','none');


        for i = 1:length(pr.stimSizeDegsDomain)
            for j = 1:length(pr.focalRegionDomain)
                for k = 1:length(pr.focalPropLListDomain)

                    [theConeMosaic, mosaicConeInfo] = setConeProportions(pr.focalRegionDomain(j), ...
                        pr.focalPropLListDomain(k), pr.focalVariantDomain(h), theConeMosaic, pr.eccXDegs, pr.eccYDegs, ...
                        pr.stimSizeDegsDomain(i), pr.fieldSizeMinutes, pr.regionVariant, pr.propL, pr.propS);

                    % Plot the mosaic in the montage
                    theAxes = nexttile;
                    figureHandle = theFig;
                    theConeMosaic.visualizeMosaic(figureHandle,theAxes);
                    set(gca, 'xticklabel', [], 'yticklabel', []);
                    set(gca, 'xlabel', [], 'ylabel', []);

                    hold on;

                    if viewBounds
                        % Pull region boundary info
                        xBounds = mosaicConeInfo.xBounds;
                        yBounds = mosaicConeInfo.yBounds;

                        % Superimpose the boundaries
                        for w = length(xBounds)
                            rectangle('Position', ...
                                [xBounds(w,1) yBounds(w,1) ...
                                (xBounds(w,2) - xBounds(w,1)) ...
                                (yBounds(w,2) - yBounds(w,1))], ...
                                'LineWidth', 3)
                        end
                    end

                    % Set plot figure information
                    if i == 1 & j == 1
                        title([num2str(pr.focalPropLListDomain(k)) 'L'])
                    end

                    if k == 1
                        ylabel({regionAxesNames(j); [num2str(pr.stimSizeDegsDomain(i)*60) ' Stim'] }, 'FontWeight', 'bold')
                    end
                end
            end
        end

        if pr.viewMosaicMontage
            set(gcf, 'Position', [595 5 1361 972]);
            title(t, [stageTitle ' Mosaic Montage, Variant ' num2str(pr.focalVariantDomain)], 'FontSize', 40)
        end
    end
end

end

