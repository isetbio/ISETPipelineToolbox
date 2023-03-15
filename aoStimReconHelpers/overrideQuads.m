function theConeMosaic = overrideQuads(theConeMosaic, eccX, eccY, chrom)
% Manually establish conditions to create a custom-formatted mosaic.
% Uses actual index position along the cone mosaic so may not generalize
% between FOV. Establish new elseif argument with new quadSeq and include
% description of mosiac before new entry.


%% Initialize global variables
global allDone; global bookKeep; 
allDone = 0;
bookKeep = cell(1,3);


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
% From Seq 9 - 29 increasing M cone count by one over the spread
% [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
%         350 459 429 349 469 416 369 351]
elseif strcmp(chrom, 'quadSeq9')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

% Quad Seq 10 - 30 arcmin
elseif strcmp(chrom, 'quadSeq10')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq11')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq12')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq13')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq14')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq15')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq16')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq17')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq18')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq19')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq20')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq21')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq22')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq23')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq24')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq25')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq26')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq27')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq28')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq29')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369 351];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()





% READ ME!! Introduced makeshift GUI callback for more direct manipulation.
% To be used only when first setting up mosaic positions.
% Before closing must capture the updated index list for each cone class and
% copy/paste those values to match the form above, then comment out the GUI
% portion. This step is necessary to avoid having to reset the cones every
% time. 
 

% Quad Seq 30 - 30 arcmin (QS 8 repeat with 3 S cones)
elseif strcmp(chrom, 'quadSeq30')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newSIndices = [404	361	370];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
    theConeMosaic.visualizeMosaic()

% Quad Seq 31 - 30 arcmin (QS 14 repeat with 3 S cones)
elseif strcmp(chrom, 'quadSeq31')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

%     theConeMosaic.visualizeMosaic(); hold on;
%     g = @(x, y) mouseClick(x, y, theConeMosaic);
%     set(gcf, 'WindowButtonDownFcn', g)
%     keyboard

    newSIndices = [404	361	370];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

% Quad Seq 32 - 30 arcmin (QS 29 repeat with 3 S cones)
elseif strcmp(chrom, 'quadSeq32')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369 351];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361	370];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()








% Quad Seq 33 - 30 arcmin (QS 29 but doubled M)
elseif strcmp(chrom, 'quadSeq33')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369 351];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    
%     theConeMosaic.visualizeMosaic(); hold on;
%     g = @(x, y) mouseClick(x, y, theConeMosaic);
%     set(gcf, 'WindowButtonDownFcn', g)
%     keyboard
    
    pad = 0.12;

    newMIndices = find(...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) <  eccX + pad & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,1) >  eccX - pad & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) <  eccY + pad & ...
        theConeMosaic.Mosaic.coneRFpositionsDegs(:,2) >  eccY - pad);
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
    theConeMosaic.visualizeMosaic()
    keyboard

% Quad Seq 34 - 30 arcmin (QS 29 but surround alternate L/S)
elseif strcmp(chrom, 'quadSeq34')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369 351];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    theConeMosaic.visualizeMosaic(); hold on;
    g = @(x, y) mouseClick(x, y, theConeMosaic);
    set(gcf, 'WindowButtonDownFcn', g)
    keyboard

    newMIndices = [newMIndices];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);   
    theConeMosaic.visualizeMosaic()
























else
    error(['Entered chrom has no alloted override sequence']);
end
end



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
