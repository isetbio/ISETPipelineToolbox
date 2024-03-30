function coneProp = calcConeProportions(pr, cnv, renderType, annWidthArc, visualizeAnnuli)
% Synopsis:
%    Function to pull out cone class proportionality under a given region
%
% Description:
%    Script to take the desired testing conditions for a given stimulation
%    and call the associated retinal mosaic (load render structure). Then
%    determines the proportion of each cone class under the stimulus as
%    well as that within a surrounding square annulus region of desired
%    width. 
%
% See also: aoStimReconRunMany, aoStimRecon

% History:
%   05/30/23  chr  Wrote as a callable function
%   03/29/24  chr  THIS FUNCTION IS NO LONGER ACTIVE, STILL HERE TO
%                  ENSURE LOSING IT DOESN'T BREAK THINGS, BUT PLAN IS TO 
%                  DELETE AT EARLIEST CONVENIENCE

%% Close existing figures
close all;

%% Grab render matrices/files

% Grab the appropriate renderStructure, in most conditions will be the same
switch renderType
    case 'forward'
        if (~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'file'))
            error('Forward render strucure not cached')
        else
            clear forwardRenderStructure;
            load(fullfile(cnv.renderDir, 'xRenderStructures', cnv.forwardRenderStructureName),'renderStructure');
            forwardRenderStructure = renderStructure; clear renderStructure;
            grabRenderStruct(forwardRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
                pr.nPixels, cnv.forwardPupilDiamMM, pr.forwardAORender, pr.forwardDefocusDiopters);
        end

        % Set forward variables from loaded/built structure. Scale matrix with display factor.
        theConeMosaic = forwardRenderStructure.theConeMosaic;
        clear forwardRenderStructure;

    case 'recon'
        % Grab recon cone mosaic and render matrix
        if (~exist(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'file'))
            error('Recon render strucure not cached');
        else
            clear reconRenderStructure;
            load(fullfile(cnv.renderDir, 'xRenderStructures', cnv.reconRenderStructureName),'renderStructure');
            reconRenderStructure = renderStructure; clear renderStructure;
            grabRenderStruct(reconRenderStructure, pr.eccXDegs, pr.eccYDegs, cnv.fieldSizeDegs, ...
                pr.nPixels, cnv.reconPupilDiamMM, pr.reconAORender, pr.reconDefocusDiopters);
        end

        % Set recon variables from loaded/built structure. Scale matrix with display factor.
        theConeMosaic = reconRenderStructure.theConeMosaic;
        clear reconRenderStructure;

end


%% Proportionality Code

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