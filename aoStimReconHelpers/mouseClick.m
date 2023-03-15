function mouseClick(source, eventdata, theConeMosaic)
% Filler fxn grabbed from MATLAB support team to display the mouse
% coordinates over a figure. To be expanded upon. 

C = get (gca, 'CurrentPoint'); 
mousePos = [C(1,1) C(1,2)];
title(gca, ['(X,Y) = (', num2str(mousePos(1)), ', ',num2str(mousePos(2)), ')']);

coneChoice = dsearchn(theConeMosaic.Mosaic.coneRFpositionsDegs, mousePos);
conePos = theConeMosaic.Mosaic.coneRFpositionsDegs(coneChoice,:);

prompt = "Enter desired cone class: ";
classChoice = input(prompt);

if classChoice == 1
    theConeMosaic.Mosaic.reassignTypeOfCones(coneChoice, cMosaic.LCONE_ID);
    newPoint = 'r.';
    plot(conePos(1), conePos(2), newPoint,'MarkerSize',50);
elseif classChoice == 2
    theConeMosaic.Mosaic.reassignTypeOfCones(coneChoice, cMosaic.MCONE_ID);
    newPoint = 'g.';
    plot(conePos(1), conePos(2), newPoint,'MarkerSize',50);
elseif classChoice == 3
    theConeMosaic.Mosaic.reassignTypeOfCones(coneChoice, cMosaic.SCONE_ID);
    newPoint = 'b.';
    plot(conePos(1), conePos(2), newPoint,'MarkerSize',50);
elseif classChoice == 0
    return
end

global allDone; global bookKeep;
allDone = input('All done? ');


if allDone
    thisFig = get(gcf, "Number");
    close(figure(thisFig))
    bookKeep(1) = {theConeMosaic.Mosaic.lConeIndices'};
    bookKeep(2) = {theConeMosaic.Mosaic.mConeIndices'};
    bookKeep(3) = {theConeMosaic.Mosaic.sConeIndices'};
    openvar('bookKeep')
    disp(["Don't forget to copy the new indices into overrideQuads and comment out the GUI portion"])
end
