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
p.addParameter('randomSLocations', false, @islogical);
p.addParameter('linearDistanceDecrement', 0.05, @isnumeric);
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
        propL = focalPropL * ones(1,3);
        regionVariant = focalVariant * ones(1,3);
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

    % Set the S cones, either randomly or using a method that ensures they
    % are spaced.%%%
    linearDistancePerSConeFactor = 1;
    runCounter = 0;

    if (p.Results.randomSLocations | ~isempty(innerCones))
        newSMosaicInd = randsample(regionCones, newSAmount)';
    else
        % We will enforce a minimum spacing between the S cones that we
        % select.  The trick is to choose a spacing that is large enough so
        % that the S cones are reasonably spaced, but not so large that it
        % cannot be achieved.  This will depend both on the size of the
        % region and the number of S cones.

        % Get region size and set minimum distance.  As part of this, apply
        % a scaling factor down from the theoretically obtained minimum
        % distnace to account for practicalities that are hard to handle
        % analytically.
        xDelta = (xBounds(i+1,2)-xBounds(i+1,1));
        yDelta = (yBounds(i+1,2)-yBounds(i+1,1));
        regionArea = xDelta*yDelta;
        regionAreaPerSCone = regionArea/newSAmount;

        while runCounter >= 0

            % The only absolute sanity check here, if hit hard lower limit
            % on distance between S cones.
            if p.Results.linearDistanceDecrement <= 0
                error('Cannot scale interval any smaller, unable to satisfy S Cones');
            end

            adjustedLinearFactor = linearDistancePerSConeFactor - (p.Results.linearDistanceDecrement * runCounter);
            linearDistancePerSCone = adjustedLinearFactor * sqrt(regionAreaPerSCone);

            % Choose S cones
            %
            % Initialize by choosing an S cone at random.
            SConePool = regionCones;
            newSMosaicInd = randsample(SConePool,1);

            % Get the remaining S cones
            for ss = 1:newSAmount-1
                % First remove from pool any cones that are not far enough away
                % from the last one we picked.  This will also get rid of the
                % one we just picked, since it's distance from itself is 0.  We
                % only need to do this removal once for each S cone we pick,
                % since all the remaining cones are guaranteed not to be too
                % close to that one.
                %
                % Start by getting location of the S cone we just picked.
                currSConeXLocation = theConeMosaic.Mosaic.coneRFpositionsDegs(newSMosaicInd(ss),1);
                currSConeYLocation = theConeMosaic.Mosaic.coneRFpositionsDegs(newSMosaicInd(ss),2);

                % Create a new pool to draw from, which is those in the current
                % pool except for those too close to the one we just picked.
                SConePoolNew = []; SConeNewIndex = 1;
                for pp = 1:length(SConePool)
                    % Get location of each cone in the pool, and distance to
                    % the S cone whose neighbors we are excluding.
                    poolSConeXLocation = theConeMosaic.Mosaic.coneRFpositionsDegs(SConePool(pp),1);
                    poolSConeYLocation = theConeMosaic.Mosaic.coneRFpositionsDegs(SConePool(pp),2);
                    checkLinearDistance = sqrt((currSConeXLocation-poolSConeXLocation)^2 + (currSConeYLocation-poolSConeYLocation)^2);

                    % If the cone is far enough away, add it to the pool we're
                    % building up.
                    if (checkLinearDistance >= linearDistancePerSCone)
                        SConePoolNew(SConeNewIndex) = SConePool(pp);
                        SConeNewIndex = SConeNewIndex + 1;
                    end
                end

                % Update the pool with the new one we just set up.
                sharedSCones = ismember(SConePool, SConePoolNew);
                SConePool = SConePool(sharedSCones);

                % Put some arbitrary integer that has no risk of being
                % called as a placeholder. Will flag below, but this method
                % gives the system a chance to course correct instead of
                % crashing with an error.
                if (isempty(SConePool) || length(SConePool) == 1)
                    SConePool = 1;
                end

                % Pick the S cone randomly from our acceptable choices, and
                % then return to the top of the loop to continue until done.
                newSMosaicInd(ss+1) = randsample(SConePool,1);
            end

            % If recognize the arbitrary flag above, should also be the
            % case that the number of selected S cones does not match the
            % desired cone, since the arbitrary index is well outside
            % center range. If so, run again, otherwise end the loop.
            targetRegionS = sum(ismember(newSMosaicInd, regionCones));
            if SConePool == 1 & targetRegionS ~= newSAmount
                runCounter = runCounter + 1;
            elseif targetRegionS == newSAmount
                runCounter = -1;
            end
        end
    end


    % Figure out which cones can be L or M
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
    mosaicConeInfo.version = 1;

    mosaicConeInfo.regionVariant(i) = regionVariant(i);
    mosaicConeInfo.regionWidths(i) = regionWidths(i) * 60;
    mosaicConeInfo.numTotal(i) = length(regionCones);

    mosaicConeInfo.numL(i) = sum(ismember(regionCones, theConeMosaic.Mosaic.lConeIndices));
    mosaicConeInfo.numM(i) = sum(ismember(regionCones, theConeMosaic.Mosaic.mConeIndices));
    mosaicConeInfo.numS(i) = sum(ismember(regionCones, theConeMosaic.Mosaic.sConeIndices));

    mosaicConeInfo.targetPropsL(i) = propL(i);
    mosaicConeInfo.targetPropsM(i) = 1 - propL(i);
    mosaicConeInfo.targetPropsS(i) = propS(i);

    mosaicConeInfo.targetPropsLOfWholeMosaic(i) = propL(i)*(1-propS(i));
    mosaicConeInfo.targetPropsMOfWholeMosaic(i) = (1 - propL(i))*(1-propS(i));
    mosaicConeInfo.targetPropsSOfWholeMosaic(i) = propS(i);

    mosaicConeInfo.achievedPropsL(i) = mosaicConeInfo.numL(i) / (mosaicConeInfo.numL(i) + mosaicConeInfo.numM(i));
    mosaicConeInfo.achievedPropsM(i) = mosaicConeInfo.numM(i) / (mosaicConeInfo.numL(i) + mosaicConeInfo.numM(i));
    mosaicConeInfo.achievedPropsS(i) = mosaicConeInfo.numS(i) / mosaicConeInfo.numTotal(i);

    mosaicConeInfo.achievedPropsLOfWholeMosaic(i) = mosaicConeInfo.numL(i) / mosaicConeInfo.numTotal(i);
    mosaicConeInfo.achievedPropsMOfWholeMosaic(i) = mosaicConeInfo.numM(i) / mosaicConeInfo.numTotal(i);
    mosaicConeInfo.achievedPropsSOfWholeMosaic(i) = mosaicConeInfo.numS(i) / mosaicConeInfo.numTotal(i);

    eval(['mosaicConeInfo.indL.region' num2str(i) ' = newLMosaicInd;']);
    eval(['mosaicConeInfo.indM.region' num2str(i) ' = newMMosaicInd;']);
    eval(['mosaicConeInfo.indS.region' num2str(i) ' = newSMosaicInd;']);
end

% Call a side function to print all the high value mosaic info in an
% accessible text file
propTxtFile(cnv, mosaicConeInfo);

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

% visualizeIndices = false;
% if visualizeIndices
%     theConeMosaic.visualizeMosaic()
%     for i=1:length(theConeMosaic.Mosaic.coneTypes)
%         txt = int2str(i);
%         t = text((theConeMosaic.Mosaic.coneRFpositionsDegs(i,1)-0.005), theConeMosaic.Mosaic.coneRFpositionsDegs(i,2),txt);
%         t.FontSize=11;
%         hold on
%     end
% end
end