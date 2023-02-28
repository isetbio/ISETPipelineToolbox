function theConeMosaic = overrideQuads(theConeMosaic, eccX, eccY, chrom)
% Manually establish conditions to create a custom-formatted mosaic.
% Uses actual index position along the cone mosaic so may not generalize
% between FOV. Establish new elseif argument with new quadSeq and include
% description of mosiac before new entry.
%

% Quad Seq 5 - 30 arcmin
if strcmp(chrom, 'quadSeq5')
    xMax = max(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    xMin = min(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    yMax = max(theConeMosaic.Mosaic.coneRFpositionsDegs(:,2));
    yMin = min(theConeMosaic.Mosaic.coneRFpositionsDegs(:,2));
    
    % Quadrant 1
    newSIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > eccX & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > eccY);
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
    newSIndices = [];
    
    newMIndices = [246 293 248 264 209 184];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newLIndices = [161 96 100 132 110];
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    
    % Quadrant 2
    newMIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > eccX - ((eccX - xMin)  / 2) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < eccX & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > eccY);
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newLIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < eccX - ((eccX - xMin)  / 2) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > xMin & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > eccY);
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    newSIndices = [632 633 620];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
    newSIndices = [];
    
    newMIndices = [750 720 763 693];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newLIndices = [512 571 514 542];
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    
    % Quadrant 3
    newMIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > eccX - ((eccX - xMin)  / 3) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < eccX & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < eccY);
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newSIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) > eccX - (2 * ((eccX - xMin)  / 3)) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < eccX - ((eccX - xMin)  / 3) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < eccY);
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
    newSIndices = [];
    
    newLIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) < eccX - (2 * ((eccX - xMin)  / 3)) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < eccY);
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    newMIndices = [726];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newLIndices = [509];
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    % Quadrant 4 background
    newLIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) <  eccX + (( eccX - xMin)  / 3) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) >  eccX & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) <  eccY & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > yMin / 2);
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    newSIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) <  eccX + (2 * (( eccX - xMin)  / 3)) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) >  eccX + (( eccX - xMin)  / 3) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) <  eccY & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > yMin / 2);
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
    newSIndices = [];
    
    newMIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) >  eccX + (2 * (( eccX - xMin)  / 3)) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) <  eccY & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) > yMin / 2);
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newLIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) <  eccX + (( eccX - xMin)  / 2) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) >  eccX & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < yMin / 2);
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];
    
    newMIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) >  eccX + (( eccX - xMin)  / 2) & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) < yMin / 2);
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newSIndices = [280 299 258 108 138 107 266 136];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
    newSIndices = [];
    
    newMIndices = [268 257 269];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    newMIndices = [];
    
    newLIndices = [127 156 105 240];
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
    newLIndices = [];

% Quad Seq 8 - 30 arcmin
elseif strcmp(chrom, 'quadSeq8')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

% Quad Seq 9 - 30 arcmin
elseif strcmp(chrom, 'quadSeq9')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);







end
end


%
% Show the mosaic with index numbers superimposed when
% creating/formatting.
% 
% theConeMosaic.visualizeMosaic()
% for i=1:length(theConeMosaic.Mosaic.coneTypes)
%     txt = int2str(i);
%     t = text((theConeMosaic.Mosaic.coneRFpositionsDegs(i,1)-0.005), theConeMosaic.Mosaic.coneRFpositionsDegs(i,2),txt);
%     t.FontSize=11;
%     hold on
% end
