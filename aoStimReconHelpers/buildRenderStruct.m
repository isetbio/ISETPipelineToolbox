function renderStructure = buildRenderStruct(pr, cnv, stage, varargin)
% Synopsis:
%    Build render matrix if desired/needed
%
% Description:
%    [EXPAND] Run this function if you would like to rebuild a new mosaic and
%    render matrix.  This also gets run if there is no cached file corresponding
%    to the desired parameters. Once built, this file can be loaded from cache
%    for quicker running.
%
% See also: aoStimRecon, aoStimReconRunMany, buildBaseMosaic,
%           buildMosaicMontage.

% History:
%   08/16/22  chr  Made it a callable function
%   08/25/22  chr  Included portion for dichromacy
%   08/26/22  dhb, chr  Convert to main file, edit cone mosaic options

%% Set the stage and build base mosaic
%
% Unpack the variables that change dependent on whether building a forward
% mosaic or recon mosaic. useDisplay in reference to whether a display
% needs to actually be loaded to achieve desired function (true for render
% struct, false for montage)
useDisplay = true;
st = unpackStage(pr, cnv, stage);
[theConeMosaic, theDisplay] = buildBaseMosaic(pr, cnv, st, useDisplay);

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
    renderStructure.fieldSizeDegs = pr.fieldSizeMinutes / 60;
    renderStructure.eccX = pr.eccXDegs;
    renderStructure.eccY = pr.eccYDegs;
    renderStructure.nPixels = pr.nPixels;
    renderStructure.pupilDiamMM = st.pupilDiamMM;
    renderStructure.AORender = st.aoRender;
    renderStructure.defocusDiopters = st.defocusDiopters;
    renderStructure.mosaicConeInfo = mosaicConeInfo;

    % Save it to the output location
    save(fullfile(st.renderDirFull, cnv.renderName),'renderStructure','-v7.3');
else
    %% Build render matrix based on edited mosaic
    %
    % Build the custom mosaics using setConeProportions.  We add the
    % expansion parameter here so that the cone proportions are set within
    % the expanded region.
    [theConeMosaic, mosaicConeInfo] = setConeProportions(pr.focalRegion, ...
        pr.focalPropL, pr.focalVariant, theConeMosaic, pr.eccXDegs, pr.eccYDegs, ...
        pr.stimSizeDegs+pr.forwardOpticalBlurStimSizeExpansionDegs, ...
        pr.fieldSizeMinutes, pr.regionVariant, pr.propL, pr.propS, ...
        'annulusWidthArc',pr.annulusWidthArc);

    % Build the render structure for the custom mosaic.
    %
    % Note that the "forwardRender" method is a generic
    % method that builds the render matrix for forward
    % computations. This use of the word "forward" is
    % separate from the distinction we make between
    % "forward" and "recon" matrices.  A little unfortunate
    % but we'll live with the potential confusion.
    theConeMosaic.Display = theDisplay;
    renderMatrix = theConeMosaic.forwardRender([pr.nPixels pr.nPixels 3], ...
        true, true, 'useDoublePrecision', true);
    renderMatrix = double(renderMatrix);

    % Push new info back into structure and save
    renderStructure.theDisplay = theDisplay;
    renderStructure.renderMatrix = renderMatrix;
    renderStructure.theConeMosaic = theConeMosaic;
    renderStructure.fieldSizeDegs = pr.fieldSizeMinutes / 60;
    renderStructure.eccX = pr.eccXDegs;
    renderStructure.eccY = pr.eccYDegs;
    renderStructure.nPixels = pr.nPixels;
    renderStructure.pupilDiamMM = st.pupilDiamMM;
    renderStructure.AORender = st.aoRender;
    renderStructure.defocusDiopters = st.defocusDiopters;
    renderStructure.mosaicConeInfo = mosaicConeInfo;

    % Save it to the output location
    save(fullfile(st.renderDirFull, cnv.renderName),'renderStructure','-v7.3');

    %% Create an excel file with all the relevant information
    %
    % Create the filename
    filename = fullfile(st.renderDirFull, ...
        [cnv.renderName, '_SummaryInfo.xlsx']);

    % Build a table with the pertinent info
    A = [...
        mosaicConeInfo.regionVariant; ...
        mosaicConeInfo.regionWidths; ...
        mosaicConeInfo.numTotal; ...

        mosaicConeInfo.numL; ...
        mosaicConeInfo.numM; ...
        mosaicConeInfo.numS; ...

        mosaicConeInfo.targetPropsL; ...
        mosaicConeInfo.targetPropsM; ...
        mosaicConeInfo.targetPropsS; ...

        mosaicConeInfo.achievedPropsL; ...
        mosaicConeInfo.achievedPropsM; ...
        mosaicConeInfo.achievedPropsS; ...

        mosaicConeInfo.targetPropsLOfWholeMosaic; ...
        mosaicConeInfo.targetPropsMOfWholeMosaic; ...
        mosaicConeInfo.targetPropsSOfWholeMosaic; ...

        mosaicConeInfo.achievedPropsLOfWholeMosaic; ...
        mosaicConeInfo.achievedPropsMOfWholeMosaic; ...
        mosaicConeInfo.achievedPropsSOfWholeMosaic; ...
        ];

    % Format the array to table
    T = array2table(A, 'VariableNames', {'Center' 'Near Surround' 'Distant Surround'}, ...
        'RowNames', {'Region Variant' 'Region Widths' 'Total Cones' ...
        'Num L' 'Num M' 'Num S' ...
        'Target Prop L' 'Target Prop M' 'Target Prop S' ...
        'Achieved Prop L' 'Achieved Prop M' 'Achieved Prop S' ...
        'Target Prop L Whole' 'Target Prop M Whole' 'Target Prop S Whole' ...
        'Achieved Prop L Whole' 'Achieved Prop M Whole' 'Achieved Prop S Whole'});

    % Save the output
    writetable(T, filename, "Sheet", 1);
end
end
