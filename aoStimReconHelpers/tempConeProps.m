% Not a function, moreso an executable script while calling a greater
% function to create mosaics. Should be done by clicking the "Run Section"
% button up above 




annWidthArc = [2; 2];
visualizeAnnuli = true;
pr.stimSizeDegs = 3.5/60;

seqNum = 142;
setProps = false;
justView = true;
setPropsAnn = false;

fixedS = true;
posS = [404 361];

imageExp = true; 


if setProps
    percentL = [0.5];
    percentS = [0.0];
end

if setPropsAnn
    annulusL = [1.0];
    annulusS = [0.0];
end

if justView
    newLIndices = storedQuadIndices{1,seqNum}{1};
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = storedQuadIndices{1,seqNum}{2};
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = storedQuadIndices{1,seqNum}{3};
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
end


% Initialize the structure to hold the information
coneProp = struct;

% For each series of annuli conditions
for j = 1:size(annWidthArc, 1)

    % Initialize helpful variables
    innerCones = [];
    regionCones = [];

    % Round the stimulus size to 2 decimal points for convenience
    spanDegs = round(pr.stimSizeDegs, 2);

    % Convert the provided field size into degrees
    fieldSizeDegs = pr.fieldSizeMinutes/60;

    % Establish brackets to note if desired annulus exceeds field size.
    xCheck = [pr.eccXDegs - fieldSizeDegs/2 pr.eccXDegs + fieldSizeDegs/2];
    yCheck = [pr.eccYDegs - fieldSizeDegs/2 pr.eccYDegs + fieldSizeDegs/2];

    % Establish the first inner bound. Default is to set the first inner bound
    % equal to the center of the FOV
    xInnerBounds = [pr.eccXDegs pr.eccXDegs];
    yInnerBounds = [pr.eccYDegs pr.eccYDegs];
    % Should include an option to override such that central region is different from
    % stimulus overlay

    % Create annulus list starting with the centermost point and extending
    % outward based on input annulus width (given in arcmin for convenience but
    % converted to degrees). Final value is the full field
    annWidthDeg = [spanDegs/2 annWidthArc(j,:)/60 fieldSizeDegs/2];

    % Include optional component to visualize the annuli borders
    if visualizeAnnuli
        theConeMosaic.visualizeMosaic()
    end

    % For each desired annulus
    for i = 1:size(annWidthDeg, 2)

        % Establish the new x and y bounds as some extension beyond the inner
        % bounds
        xBounds = [xInnerBounds(1)-annWidthDeg(i), xInnerBounds(2)+annWidthDeg(i)];
        yBounds = [yInnerBounds(1)-annWidthDeg(i), yInnerBounds(2)+annWidthDeg(i)];

        % Superimpose image of rectangle on retinal mosaic defining annulus
        if visualizeAnnuli
            rectangle('Position', [xBounds(1) yBounds(1) (xBounds(2) - xBounds(1)) (yBounds(2) - yBounds(1))], 'LineWidth', 3)
        end

        % Find all the cones centered in the bounded square region
        regionCones = find(...
            theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > xBounds(1) & ...
            theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > yBounds(1) & ...
            theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < xBounds(2) & ...
            theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < yBounds(2));

        % Select for cones that do not also occupy any inner/previously
        % calculated annuli. Not applied to the final annWidth value to
        % allow capture of full FOV mosaic proportions
        if i ~= size(annWidthDeg, 2)
            regionCones = regionCones(~ismember(regionCones, innerCones));
        end




        % THIS IS THE PART OF THE CODE THAT DOES THE WHOLE CHANGE TO THIS PERCENTAGE INSTEAD.

        if setProps
            if i == 1
                % Initialize a index tracker using boolean vectors for
                % efficiency in L/M random assignments
                indTracker = false(1,length(regionCones));

                % Convert the desired L percentage to a number of cones
                newLAmount = round(length(regionCones) * percentL);

                % Select the desired number of cones from the region to be
                % converted to L cones, represented as true in the tracker
                newLRegionInd = randperm(length(regionCones), newLAmount);
                indTracker(newLRegionInd) = true;

                % Convert the desired S percentage to a number of cones and
                % apply this number to the mosaic regions
                newSAmount = round(length(regionCones) * percentS);
                newSRegionInd = randperm(length(regionCones), newSAmount);

                % Apply desired percentages to the L and S cones using the
                % tracker. Then override all with aplied S cone values
                newLMosaicInd = regionCones(indTracker);
                newMMosaicInd = regionCones(~indTracker);
                newSMosaicInd = regionCones(newSRegionInd);

                % Apply if want the S cones to stay fixed in place.
                if fixedS
                    newSMosaicInd = posS;
                end

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





            end
        end











        if setPropsAnn
            if i == 2
                % Initialize a index tracker using boolean vectors for
                % efficiency in L/M random assignments
                indTracker = false(1,length(regionCones));

                % Convert the desired L percentage to a number of cones
                newLAmount = round(length(regionCones) * annulusL);

                % Select the desired number of cones from the region to be
                % converted to L cones, represented as true in the tracker
                newLRegionInd = randperm(length(regionCones), newLAmount);
                indTracker(newLRegionInd) = true;

                % Convert the desired S percentage to a number of cones and
                % apply this number to the mosaic regions
                newSAmount = round(length(regionCones) * annulusS);
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





            end
        end








        % Determine the total number of cones in the region, as well as
        % class-specific numbers.
        numAll = length(regionCones);
        numL = sum(ismember(regionCones, theConeMosaic.Mosaic.lConeIndices));
        numM = sum(ismember(regionCones, theConeMosaic.Mosaic.mConeIndices));
        numS = sum(ismember(regionCones, theConeMosaic.Mosaic.sConeIndices));

        % Store values in an appropriate struct. The following format
        % allows names to be iteratively changed within the loop
        eval(['coneProp.series' num2str(j) '.size( ' num2str(i) ') = annWidthDeg(' num2str(i) ') * 60;']);
        eval(['coneProp.series' num2str(j) '.numAll( ' num2str(i) ') = numAll;']);
        eval(['coneProp.series' num2str(j) '.numL( ' num2str(i) ') = numL;']);
        eval(['coneProp.series' num2str(j) '.numM( ' num2str(i) ') = numM;']);
        eval(['coneProp.series' num2str(j) '.numS( ' num2str(i) ') = numS;']);

        % Update such that annuli border become the bounds and cones now become
        % inner limits for the next annulus layer
        if i ~= (size(annWidthDeg, 2)-1)
            xInnerBounds = xBounds;
            yInnerBounds = yBounds;
            innerCones = [innerCones; regionCones];

            % Else portion is because going into the final annWidthDeg value we
            % don't want an inner bound, instead want to capture full FOV
        else
            xInnerBounds = [pr.eccXDegs pr.eccXDegs];
            yInnerBounds = [pr.eccYDegs pr.eccYDegs];
        end

    end
end



coneProp.series1
coneProp.series2


if setProps | setPropsAnn
    bookKeep(1) = {theConeMosaic.Mosaic.lConeIndices'};
    bookKeep(2) = {theConeMosaic.Mosaic.mConeIndices'};
    bookKeep(3) = {theConeMosaic.Mosaic.sConeIndices'};
    storedQuadIndices(1,(seqNum)) = {bookKeep};
end


if imageExp
    f = gcf;
    exportgraphics(f, [num2str(seqNum) '.png'], 'Resolution', 600);
end

% save(fullfile(storedDir, "storedQuadIndices.mat"), "storedQuadIndices")
