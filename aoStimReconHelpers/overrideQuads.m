function theConeMosaic = overrideQuads(theConeMosaic, eccX, eccY)
%% QuadSeq5

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


% % Show it w/ numbers
% if (indDisp)
%     theConeMosaic.visualizeMosaic()
%     for i=1:length(inputMosaic.coneTypes)
%         txt = int2str(i);
%         t = text((theConeMosaic.coneRFpositionsDegs(i,1)-0.005), theConeMosaic.Mosaic.coneRFpositionsDegs(i,2),txt);
%         t.FontSize=11;
%         hold on
%     end
% end
