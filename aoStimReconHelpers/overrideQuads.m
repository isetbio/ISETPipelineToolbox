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
% time. Prone to crashing so if switching a lot of cones do so in chunks 
 

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

    newSIndices = [404	361	370];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

%     theConeMosaic.visualizeMosaic(); hold on;
%     g = @(x, y) mouseClick(x, y, theConeMosaic);
%     set(gcf, 'WindowButtonDownFcn', g)
%     keyboard

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
    

% Quad Seq 34 - 30 arcmin (QS 29 but surround alternate L/S)
elseif strcmp(chrom, 'quadSeq34')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369 351];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [3	4	7	8	10	11	16	25	33	39	42	43	45	48	53	54	61	64	65	67	69	74	80	86	87	88	90	95 ...
    	102	104	106	107	109	111	113	117	119	122	129	132	142	151	152	154	155	157	165	173	174	178	179	181	183	184	188	193	194	200	203	...
        210	212	213	219	221	223	225	231	242	247	254	256	257	260	263	266	270	272	275	276	277	279	281	293	295	296	298	299	305	306	308	...
        310	313	320	321	325	331	337	339	343	345	348	352	356	363	364	388	397	398	399	402	405	406	407	409	418	425	427	433	436	440	442	...
        444	460	461	463	464	473	475	477	489	490	491	493	495	497	498	506	510	520	524	526	528	531	539	544	551	554	557	560	564	571	573	...
        575	579	589	592	597	600	602	605	610	611	624	629	631	633	640	644	646	647	654	656	660	666	671	673	674	675	680	695	697	698	699	...
        712	717	718	720	722	727	728	730	731	733	736	739	748	749	750	751	752	754	759	767	770	772	777	780	782	785	789	791	792	794	795	...
        807	814	817	825	827	828	830	831	834	836	357	358	160	435	534	487	548	616	648	756	826];

    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

    
%     theConeMosaic.visualizeMosaic(); hold on;
%     g = @(x, y) mouseClick(x, y, theConeMosaic);
%     set(gcf, 'WindowButtonDownFcn', g)
%     keyboard

    







% Quad Seq 8 - 30 arcmin
elseif strcmp(chrom, 'quadSeq35')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()


% Quad Seq 10 - 30 arcmin
elseif strcmp(chrom, 'quadSeq36')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq37')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq38')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq39')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()
    
elseif strcmp(chrom, 'quadSeq40')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()
    
elseif strcmp(chrom, 'quadSeq41')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()
    
elseif strcmp(chrom, 'quadSeq42')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()

elseif strcmp(chrom, 'quadSeq43')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()
    
elseif strcmp(chrom, 'quadSeq44')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
    theConeMosaic.visualizeMosaic()
    
elseif strcmp(chrom, 'quadSeq45')
    newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
    theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

    newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        350 459 429 349 469 416 369 351];
    theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

    newSIndices = [404	361];
    theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);   
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
