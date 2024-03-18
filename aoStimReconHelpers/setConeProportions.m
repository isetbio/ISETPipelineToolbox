function [theConeMosaic, mosaicConeInfo] = setConeProportions(focalRegion, ...
    focalPropL, focalVariant, theConeMosaic, eccXDegs, eccYDegs, ...
    stimSizeDegs, fieldSizeMinutes, regionVariant, propL, propS, varargin)
% Assign cone proportions and permutation variants for each region
%
% Description:
%    When building a new coneMosaic, use this function to apply desired
%    L cone proportionality to a focal region. Options include: center,
%    nearSurround, distantSurround, multiple, and global (the full mosaic).
%    focalVariant alters the cone positions within a given proportionality
%    and region.
%
% See also: buildRenderStruct.m, aoStimReconRunMany_small_quads_chr.m
%
% History:
%   03/07/24  chr  Make it a callable function

p = inputParser;
p.addParameter('viewMosaic', false, @islogical);
p.addParameter('annulusWidthArc', 2, @isnumeric);
p.addParameter('stimCenter', true, @islogical)
p.addParameter('centerSizeDegs', [], @isnumeric);
parse(p, varargin{:});

%% Initialize variables
% Initialize the structure to hold the information
mosaicConeInfo = struct;

% Initialize helpful variables
innerCones = [];
regionCones = [];
xBounds = [];
yBounds = [];

% Establish the innermost bound as the FOV center
xBounds(1,:) = [eccXDegs eccXDegs];
yBounds(1,:) = [eccYDegs eccYDegs];

%% Adjust focal region
% Override default values with desired region proportion and variant 
switch focalRegion
    case 'center'
        propL(1) = focalPropL;
        regionVariant(1) = focalVariant;
    case 'nearSurround'
        propL(2) = focalPropL;
        regionVariant(2) = focalVariant;
    case 'distantSurround'
        propL(3) = focalPropL;
        regionVariant(3) = focalVariant;
    case 'multiple'
        if length(focalPropL) ~= 3 || length(focalVariant) ~=3
            error(['Must provide proportions and variants for all three regions'])
        end
        propL = focalPropL;
        regionVariant = focalVariant;
    case 'global'
        propL = focalPropL;
        regionVariant = focalVariant;
    otherwise
        error(['Unrecognized focal region entered'])
end

%% Set the boundaries
% Establish boundaries for each focal region unless making global changes
switch focalRegion
    case 'global'
        fieldSizeDegs = fieldSizeMinutes/60;
        regionWidths = fieldSizeDegs/2;
    otherwise
        if p.Results.stimCenter
            centerWidth = stimSizeDegs/2;
        else
            centerWidth = p.Results.centerSizeDegs/2;
        end

        nearSurroundWidth = p.Results.annulusWidthArc/60;
        fieldSizeDegs = fieldSizeMinutes/60;
        distantSurroundWidth = (fieldSizeDegs/2) - centerWidth - nearSurroundWidth;

        % Establish a list of values corresponding to each of the regions
        regionWidths = [centerWidth nearSurroundWidth distantSurroundWidth];
end

%% Cycle through the cone assignments for each region
% For each desired region
for i = 1:length(regionWidths)

    % Expand outward by adding the next width to the current border
    xBounds(i+1,:) = [xBounds(i,1)-regionWidths(i), xBounds(i,2)+regionWidths(i)];
    yBounds(i+1,:) = [yBounds(i,1)-regionWidths(i), yBounds(i,2)+regionWidths(i)];

    % Find all the cones centered in the bounded square region
    regionCones = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > xBounds(i+1,1) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > yBounds(i+1,1) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < xBounds(i+1,2) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < yBounds(i+1,2));

    % Select for cones that do not also occupy a previous region
    regionCones = regionCones(~ismember(regionCones, innerCones));

    % Convert the desired S percentage to a number of cones and
    % apply this number to the mosaic regions. Assign the random generator
    % to ensure consistency and ID unused remaining cones
    randState = rng(regionVariant(i));
    newSAmount = round(length(regionCones) * propS(i));
    newSMosaicInd = randsample(regionCones, newSAmount)';
    conesRemaining = regionCones(~ismember(regionCones, newSMosaicInd));

    % With the remaining cones, assign the desired L percentage, again
    % assigning the rng. Any cones left unselected will be converted to M
    randState = rng(regionVariant(i));
    newLAmount = round(length(conesRemaining) * propL(i));
    newLMosaicInd = randsample(conesRemaining, newLAmount)';
    newMMosaicInd = conesRemaining(~ismember(conesRemaining, newLMosaicInd))';

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

    % Update running list of which cones have been edited
    innerCones = [innerCones; regionCones];

    % Store relevant info in a structure
    mosaicConeInfo.targetPropsL(i) = propL(i);
    mosaicConeInfo.targetPropsM(i) = 1 - propL(i) - propS(i);
    mosaicConeInfo.targetPropsS(i) = propS(i);

    mosaicConeInfo.regionVariant(i) = 1; %%%%%%%%%%;
    mosaicConeInfo.regionWidths(i) = regionWidths(i) * 60;
    mosaicConeInfo.numTotal(i) = length(regionCones);

    mosaicConeInfo.numL(i) = sum(ismember(regionCones, theConeMosaic.Mosaic.lConeIndices));
    mosaicConeInfo.numM(i) = sum(ismember(regionCones, theConeMosaic.Mosaic.mConeIndices));
    mosaicConeInfo.numS(i) = sum(ismember(regionCones, theConeMosaic.Mosaic.sConeIndices));

    mosaicConeInfo.achievedPropsL(i) = mosaicConeInfo.numL(i) / mosaicConeInfo.numTotal(i);
    mosaicConeInfo.achievedPropsM(i) = mosaicConeInfo.numM(i) / mosaicConeInfo.numTotal(i);
    mosaicConeInfo.achievedPropsS(i) = mosaicConeInfo.numS(i) / mosaicConeInfo.numTotal(i);

    eval(['mosaicConeInfo.indL.region' num2str(i) ' = newLMosaicInd;']);
    eval(['mosaicConeInfo.indM.region' num2str(i) ' = newMMosaicInd;']);
    eval(['mosaicConeInfo.indS.region' num2str(i) ' = newSMosaicInd;']);
end

% Save the boundary information as well
mosaicConeInfo.xBounds = xBounds; 
mosaicConeInfo.yBounds = yBounds;

%% Visualize mosaic
% View the created mosaic with regions superimposed
if p.Results.viewMosaic
    theConeMosaic.visualizeMosaic(); hold on
    for w = 2:length(xBounds)
        rectangle('Position', [xBounds(w,1) yBounds(w,1) (xBounds(w,2) - xBounds(w,1)) (yBounds(w,2) - yBounds(w,1))], 'LineWidth', 3)
    end
end
end