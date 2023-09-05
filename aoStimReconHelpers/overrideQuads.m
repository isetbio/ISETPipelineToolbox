function theConeMosaic = overrideQuads(theConeMosaic, eccX, eccY, chrom, storedDir, pr, cnv)
% Manually establish conditions to create a custom-formatted mosaic.
% Uses actual index position along the cone mosaic so may not generalize
% between FOV. Establish new elseif argument with new quadSeq and include
% description of mosiac before new entry.


%% Initialize global variables
global allDone; global bookKeep;
allDone = 0;
bookKeep = cell(1,3);
load(fullfile(storedDir, "storedQuadIndices.mat"));
chrom = convertStringsToChars(chrom);

% Quad Seq 5 - 30 arcmin
switch chrom
    case 'quadSeq5'
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
    case 'quadSeq8'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);




        % Quad Seq 9 - 30 arcmin
        % From Seq 9 - 29 increasing M cone count by one over the spread
        % [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
        %         350 459 429 349 469 416 369 351]
    case 'quadSeq9'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        % Quad Seq 10 - 30 arcmin
    case 'quadSeq10'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq11'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq12'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq13'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq14'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq15'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq16'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq17'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq18'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq19'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq20'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq21'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq22'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq23'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq24'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq25'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq26'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq27'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq28'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq29'
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
    case 'quadSeq30'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newSIndices = [404	361	370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

        % Quad Seq 31 - 30 arcmin (QS 14 repeat with 3 S cones)
    case 'quadSeq31'
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
    case 'quadSeq32'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361	370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

        % Quad Seq 33 - 30 arcmin (QS 29 but doubled M)
    case 'quadSeq33'
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
    case 'quadSeq34'
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
    case 'quadSeq35'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()


        % Quad Seq 10 - 30 arcmin
    case 'quadSeq36'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq37'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq38'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq39'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq40'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq41'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq42'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq43'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq44'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq45'
        newLIndices = find(theConeMosaic.Mosaic.coneRFpositionsDegs(:,1));
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()



        % Quad Seq 8 - 30 arcmin
    case 'quadSeq46'

        newLIndices = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303	304	305	306	307	308	309	310	311	312	313	314	315	316	317	318	319	320	321	322	323	324	325	326	327	328	329	330	331	332	333	334	335	336	337	338	339	340	341	342	343	344	345	346	347	348	349	350	351	352	353	354	355	356	357	358	359	360	361	362	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377	378	379	380	381	382	383	384	385	386	387	388	389	390	391	392	393	394	395	396	397	398	399	400	401	402	403	404	405	406	407	408	409	410	411	412	413	414	415	416	417	418	419	420	421	422	423	424	425	426	427	428	429	430	431	432	433	434	435	436	437	438	439	440	441	442	443	444	445	446	447	448	449	450	451	452	453	454	455	456	457	458	459	460	461	462	463	464	465	466	467	468	469	470	471	472	473	474	475	476	477	478	479	480	481	482	483	484	485	486	487	488	489	490	491	492	493	494	495	496	497	498	499	500	501	502	503	504	505	506	507	508	509	510	511	512	513	514	515	516	517	518	519	520	521	522	523	524	525	526	527	528	529	530	531	532	533	534	535	536	537	538	539	540	541	542	543	544	545	546	547	548	549	550	551	552	553	554	555	556	557	558	559	560	561	562	563	564	565	566	567	568	569	570	571	572	573	574	575	576	577	578	579	580	581	582	583	584	585	586	587	588	589	590	591	592	593	594	595	596	597	598	599	600	601	602	603	604	605	606	607	608	609	610	611	612	613	614	615	616	617	618	619	620	621	622	623	624	625	626	627	628	629	630	631	632	633	634	635	636	637	638	639	640	641	642	643	644	645	646	647	648	649	650	651	652	653	654	655	656	657	658	659	660	661	662	663	664	665	666	667	668	669	670	671	672	673	674	675	676	677	678	679	680	681	682	683	684	685	686	687	688	689	690	691	692	693	694	695	696	697	698	699	700	701	702	703	704	705	706	707	708	709	710	711	712	713	714	715	716	717	718	719	720	721	722	723	724	725	726	727	728	729	730	731	732	733	734	735	736	737	738	739	740	741	742	743	744	745	746	747	748	749	750	751	752	753	754	755	756	757	758	759	760	761	762	763	764	765	766	767	768	769	770	771	772	773	774	775	776	777	778	779	780	781	782	783	784	785	786	787	788	789	790	791	792	793	794	795	796	797	798	799	800	801	802	803	804	805	806	807	808	809	810	811	812	813	814	815	816	817	818	819	820	821	822	823	824	825	826	827	828	829	830	831	832	833	834	835	836	837	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newSIndices = [110	698	670	767	115	421	337	145	478	503	178	430	818	405	573	339	29	241	658	284	69	418	300	603	428	655	664	154	101	665	516	14	723	415	438	487	610	686	307	68	586	265	669	296	820	140	103	696	35	542	579	345	299	770	313	787	123	255	245	697	193	242	318	549	112	674	65	356	24	580	538	165	521	427	650	765	688	320	273	372	195	704	353	192];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()


        % Quad Seq 10 - 30 arcmin
    case 'quadSeq47'
        newLIndices = [1	2	3	4	5	6	7	8	10	11	12	13	14	15	16	17	18	19	20	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	39	40	41	43	44	45	46	47	48	50	51	52	53	54	55	56	57	58	59	60	61	62	63	65	66	67	69	70	71	72	73	74	75	77	78	79	80	81	82	83	84	85	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	114	115	118	119	120	121	122	123	124	126	127	128	131	133	134	135	136	137	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	156	157	158	159	160	161	162	163	164	165	166	167	168	170	171	172	173	174	175	176	177	178	179	181	182	183	184	185	186	187	188	189	190	191	193	194	195	196	197	198	200	201	202	204	205	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	241	243	244	245	246	248	249	250	251	252	253	254	256	257	258	259	260	261	262	263	264	265	266	268	269	270	271	272	273	274	276	278	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303	304	305	306	307	308	309	310	311	313	314	315	317	318	319	320	321	322	323	326	327	330	331	332	334	335	336	337	338	339	340	341	342	343	344	345	346	347	348	349	350	351	352	353	354	355	356	357	358	359	360	361	362	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377	379	380	381	382	384	386	387	389	390	391	392	393	394	395	396	397	398	399	400	401	402	403	404	405	407	408	409	411	412	413	414	415	416	417	418	419	420	421	422	423	424	425	426	427	428	429	430	431	432	433	434	435	437	438	439	440	441	442	443	445	446	447	449	450	451	452	453	454	455	456	457	458	459	460	461	462	463	464	465	466	467	468	469	470	471	473	474	475	476	477	478	479	480	481	482	483	484	485	486	487	488	489	490	491	492	493	494	495	496	497	499	500	501	502	503	504	505	506	507	508	509	510	511	512	514	515	516	517	518	519	520	521	522	523	524	525	526	527	528	529	530	532	533	534	535	536	537	538	539	540	541	542	543	544	545	546	547	548	549	550	551	552	553	555	556	557	558	559	560	561	562	563	564	565	567	568	569	570	571	572	573	574	575	576	577	578	579	580	581	582	583	584	585	586	587	588	589	590	591	592	593	594	595	596	597	598	599	600	602	603	604	605	606	607	608	609	610	611	612	613	614	615	616	617	618	619	620	621	623	624	625	627	628	629	630	631	632	633	634	635	636	637	639	641	642	644	645	646	647	648	649	650	651	652	653	654	655	657	658	659	660	661	662	663	664	665	666	668	669	670	671	672	674	675	676	677	678	679	680	681	682	683	684	685	686	687	688	689	690	691	693	694	695	696	697	698	699	700	701	702	703	704	705	706	707	708	709	710	711	713	714	715	716	717	719	720	721	722	723	724	726	727	728	729	730	731	732	733	734	735	736	737	738	739	740	741	742	744	745	746	747	748	750	751	752	754	755	756	757	758	759	760	761	762	763	764	765	766	767	768	769	771	773	774	776	777	778	779	780	781	782	784	785	786	787	788	790	791	792	793	794	795	796	797	798	799	801	803	804	805	806	807	808	809	810	811	812	813	814	815	816	817	818	819	821	822	823	824	825	826	828	829	830	832	833	834	835	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [9	21	38	42	49	64	68	76	86	112	113	116	117	125	129	130	132	138	139	155	169	180	192	199	203	206	239	240	242	247	255	267	275	277	279	312	316	324	325	328	329	333	378	383	385	388	406	410	436	444	448	472	498	513	531	554	566	601	622	626	638	640	643	656	667	673	692	712	718	725	743	749	753	770	772	775	783	789	800	802	820	827	831	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [383	41	618	32	796	619	780	427	200	216	629	822	295	622	91	492	354	600	215	78	369	523	108	370	531	673	251	327	717	567	196	614	235	224	5	301	351	244	233	194	748	686	317	382	449	389	214	783	146	680	26	262	589	506	133	746	425	197	452	714	697	375	344	242	43	583	102	275	305	681	17	648	221	192	374	557	155	165	743	451	231	733	678	144];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq48'
        newLIndices = [1	2	3	4	5	7	8	9	11	12	13	14	15	16	17	18	19	20	21	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	40	41	42	45	46	47	48	49	50	53	54	55	56	57	58	59	60	61	62	63	64	65	66	69	70	71	73	74	75	77	78	79	81	84	85	86	87	89	90	91	92	93	96	97	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	118	119	120	121	122	126	127	131	132	133	134	135	136	137	139	140	141	145	148	149	150	151	152	156	157	158	159	160	161	162	163	164	165	166	167	168	170	171	174	176	177	178	180	181	182	183	184	185	186	187	188	190	191	192	193	194	195	196	197	198	200	203	204	205	206	207	208	209	211	212	213	214	217	218	219	220	221	222	224	225	226	229	230	232	234	235	236	238	239	240	241	242	243	244	246	247	249	251	252	253	254	255	256	257	258	259	260	261	263	264	265	267	268	269	270	273	276	277	278	280	282	283	284	285	286	287	289	291	292	293	294	295	296	297	298	300	301	302	304	305	307	308	309	310	311	313	315	317	318	319	320	321	322	323	325	326	327	328	329	330	331	332	334	335	336	337	338	339	340	341	342	343	345	346	347	348	350	351	352	355	356	357	360	362	363	364	365	366	367	370	371	375	376	377	379	380	381	382	383	384	385	386	387	388	389	390	391	392	394	395	396	397	398	399	400	401	402	403	404	405	406	407	408	409	410	411	412	413	414	415	416	417	418	419	421	423	424	425	427	428	429	430	432	434	435	437	438	440	441	442	444	445	446	447	448	449	451	452	453	454	455	457	459	460	461	464	465	466	467	468	470	471	472	473	475	477	478	479	480	481	482	483	484	485	486	488	489	490	491	492	494	495	496	497	498	499	500	502	503	504	506	507	508	509	510	511	512	514	515	516	517	518	519	520	521	522	523	524	525	526	528	530	532	534	535	536	537	538	539	540	541	542	543	544	545	546	548	549	550	551	553	554	555	556	557	559	560	561	564	565	566	568	569	570	572	573	574	575	576	577	578	579	581	582	583	584	585	586	587	588	590	591	592	593	594	595	596	597	598	600	602	603	604	605	606	607	608	609	611	613	615	616	617	618	619	620	622	623	625	626	627	629	630	631	632	633	635	636	637	638	639	640	642	643	644	645	646	647	648	649	650	651	652	654	655	656	657	658	659	660	661	662	663	664	666	667	668	669	671	672	673	674	675	677	679	681	683	684	685	687	688	689	690	691	692	693	694	695	696	697	698	699	700	701	702	703	706	707	708	711	712	713	714	715	716	717	718	720	721	722	725	727	728	730	731	732	733	734	735	736	737	739	741	742	743	745	746	748	749	750	752	753	754	755	756	758	759	760	761	762	764	765	766	767	768	769	770	771	772	773	774	775	776	777	778	779	780	781	783	785	786	787	788	789	790	791	792	793	795	796	797	798	799	800	802	803	804	806	808	809	810	811	813	814	815	816	817	818	821	822	823	824	826	827	828	829	830	831	832	833	834	835	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [6	10	22	39	43	44	51	52	67	68	72	76	80	82	83	88	94	95	98	117	123	124	125	128	129	130	138	142	143	144	146	147	153	154	155	169	172	173	175	179	189	199	201	202	210	215	216	223	227	228	231	233	237	245	248	250	262	266	271	272	274	275	279	281	288	290	299	303	306	312	314	316	324	333	344	349	353	354	358	359	361	368	369	372	373	374	378	393	420	422	426	431	433	436	439	443	450	456	458	462	463	469	474	476	487	493	501	505	513	527	529	531	533	547	552	558	562	563	567	571	580	589	599	601	610	612	614	621	624	628	634	641	653	665	670	676	678	680	682	686	704	705	709	710	719	723	724	726	729	738	740	744	747	751	757	763	782	784	794	801	805	807	812	819	820	825	836	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [559	15	101	796	814	26	411	717	202	692	674	521	2	314	746	560	312	519	200	468	804	695	232	557	292	803	69	203	658	812	430	647	596	115	353	282	384	471	117	724	511	130	451	741	622	544	370	206	450	197	252	831	696	624	726	141	405	490	713	518	303	575	635	466	66	780	42	407	92	293	625	188	678	546	289	190	193	584	38	827	801	566	739	290];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq49'
        newLIndices = [2	3	4	5	6	8	9	10	12	13	14	15	17	18	19	20	21	22	23	26	28	29	30	32	33	34	36	37	38	39	40	41	42	43	44	47	48	49	52	53	54	55	56	57	61	62	63	64	65	66	67	68	69	70	71	72	73	74	77	78	80	83	84	85	87	89	90	92	96	97	98	99	101	103	104	105	106	110	111	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	130	131	133	134	135	136	137	142	144	149	150	151	152	153	154	156	159	160	161	166	169	170	171	172	173	178	179	180	181	182	183	184	185	186	187	188	189	190	192	193	196	198	199	201	203	204	205	206	208	209	210	211	212	215	216	217	218	220	222	223	224	225	227	230	231	232	233	234	235	236	238	239	240	242	246	247	248	249	250	251	254	255	256	260	261	263	265	267	268	270	272	273	274	276	277	278	280	281	283	285	286	287	288	289	290	291	292	293	294	295	297	298	299	301	302	303	304	308	312	313	315	317	320	321	322	323	325	326	328	331	332	333	335	336	337	339	340	342	343	344	347	348	350	351	352	353	354	357	360	363	364	365	367	368	369	370	372	373	374	375	376	377	378	379	381	382	383	385	386	387	388	389	390	391	393	394	395	396	398	399	400	404	405	406	409	411	412	413	415	416	417	421	422	427	428	429	432	433	434	435	436	437	438	439	440	441	442	443	444	445	447	448	449	450	451	452	453	454	455	456	457	458	460	461	463	464	465	466	467	468	469	470	471	472	473	474	476	478	479	480	483	484	485	486	489	492	493	495	496	498	499	500	502	503	504	505	506	507	509	510	511	512	513	515	518	519	520	524	525	526	527	528	530	531	532	533	535	537	538	539	540	541	543	544	545	546	547	549	551	552	553	554	557	558	559	561	562	564	566	569	570	571	574	575	576	578	579	580	581	583	584	585	586	587	588	589	590	591	592	593	594	595	597	599	601	604	605	606	607	608	609	610	611	612	613	615	617	618	620	621	622	623	625	626	627	628	629	631	632	633	637	638	639	641	642	644	646	647	649	650	651	652	653	654	657	658	660	661	662	663	664	666	668	669	670	671	672	674	675	676	677	680	682	683	684	685	686	687	688	689	691	693	695	696	697	698	699	700	702	703	705	707	708	711	712	713	714	716	718	719	720	721	722	724	727	728	729	730	731	732	733	734	735	736	737	739	740	741	742	743	745	746	748	749	750	751	753	754	755	756	758	759	760	762	763	766	768	770	773	774	775	777	779	780	781	783	784	785	786	787	788	789	790	791	792	794	795	797	801	802	803	808	809	810	811	812	813	814	815	818	819	821	825	828	829	831	832	833	834	835	836	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	7	11	16	24	25	27	31	35	45	46	50	51	58	59	60	75	76	79	81	82	86	88	91	93	94	95	100	102	107	108	109	112	129	132	138	139	140	141	143	145	146	147	148	155	157	158	162	163	164	165	167	168	174	175	176	177	191	194	195	197	200	202	207	213	214	219	221	226	228	229	237	241	243	244	245	252	253	257	258	259	262	264	266	269	271	275	279	282	284	296	300	305	306	307	309	310	311	314	316	318	319	324	327	329	330	334	338	341	345	346	349	355	356	358	359	361	362	366	371	380	384	392	397	401	402	403	407	408	410	414	418	419	420	423	424	425	426	430	431	446	459	462	475	477	481	482	487	488	490	491	494	497	501	508	514	516	517	521	522	523	529	534	536	542	548	550	555	556	560	563	565	567	568	572	573	577	582	596	598	600	602	603	614	616	619	624	630	634	635	636	640	643	645	648	655	656	659	665	667	673	678	679	681	690	692	694	701	704	706	709	710	715	717	723	725	726	738	744	747	752	757	761	764	765	767	769	771	772	776	778	782	793	796	798	799	800	804	805	806	807	816	817	820	822	823	824	826	827	830	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [112	457	693	699	696	170	453	728	101	710	746	181	64	392	689	387	341	413	103	109	713	493	217	705	48	373	587	275	325	427	723	629	56	225	306	695	337	193	479	384	718	745	652	564	591	714	52	266	4	654	400	289	179	420	227	54	67	783	320	97	345	835	811	94	441	677	270	33	110	59	569	351	512	535	437	480	669	505	666	355	108	792	540	233];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq50'
        newLIndices = [2	4	6	7	8	11	12	14	17	18	19	20	22	23	24	25	26	27	28	32	34	35	36	38	39	40	42	43	45	46	47	48	49	50	51	55	56	58	62	63	65	67	68	69	74	75	76	77	78	79	80	81	83	84	85	86	87	88	92	93	95	99	101	102	105	107	108	110	115	116	117	118	120	122	123	124	125	129	130	133	135	136	137	138	139	140	141	143	144	146	147	148	149	150	151	153	154	157	159	160	161	162	168	170	176	177	178	179	180	181	183	187	188	189	195	199	200	201	202	203	209	210	211	212	213	214	215	216	217	218	220	221	222	225	226	230	233	234	236	238	239	240	242	244	245	246	247	248	251	252	253	254	256	258	259	260	261	264	268	269	270	271	273	274	275	278	279	281	283	288	289	291	292	293	294	298	299	300	305	306	309	312	314	315	317	319	320	321	323	324	326	328	329	332	335	336	337	338	339	340	342	343	344	345	347	349	350	351	353	354	355	356	362	367	368	370	373	376	377	378	379	381	382	384	388	389	390	392	393	394	396	397	400	401	402	405	406	408	409	410	411	412	415	419	422	423	424	426	427	428	429	431	432	433	434	436	437	438	439	441	442	444	446	447	448	449	450	451	452	455	456	457	458	460	461	462	467	468	469	473	476	478	479	481	483	484	490	492	498	499	500	504	506	507	508	510	511	512	513	514	515	516	517	518	519	521	522	524	525	526	527	528	529	530	532	533	534	536	537	539	540	541	542	543	545	547	548	549	550	551	553	556	558	560	561	565	566	568	570	574	578	579	582	583	585	586	587	589	590	592	593	594	596	599	600	601	602	603	605	608	610	611	615	616	617	618	619	622	623	625	626	628	630	631	632	633	634	636	637	638	639	640	643	645	646	648	649	653	654	656	658	659	661	663	666	667	668	671	672	673	675	677	678	679	681	682	684	685	686	687	688	690	691	692	693	694	695	698	701	703	706	707	708	709	710	711	712	713	714	715	717	719	720	722	723	724	726	729	730	731	732	733	736	738	739	744	745	746	749	750	752	755	756	758	759	760	761	762	763	766	767	769	770	771	772	773	775	778	779	780	781	782	784	785	786	787	791	794	795	797	798	800	802	803	804	807	810	812	813	814	815	817	818	821	822	825	827	828	832	833	834	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	3	5	9	10	13	15	16	21	29	30	31	33	37	41	44	52	53	54	57	59	60	61	64	66	70	71	72	73	82	89	90	91	94	96	97	98	100	103	104	106	109	111	112	113	114	119	121	126	127	128	131	132	134	142	145	152	155	156	158	163	164	165	166	167	169	171	172	173	174	175	182	184	185	186	190	191	192	193	194	196	197	198	204	205	206	207	208	219	223	224	227	228	229	231	232	235	237	241	243	249	250	255	257	262	263	265	266	267	272	276	277	280	282	284	285	286	287	290	295	296	297	301	302	303	304	307	308	310	311	313	316	318	322	325	327	330	331	333	334	341	346	348	352	357	358	359	360	361	363	364	365	366	369	371	372	374	375	380	383	385	386	387	391	395	398	399	403	404	407	413	414	416	417	418	420	421	425	430	435	440	443	445	453	454	459	463	464	465	466	470	471	472	474	475	477	480	482	485	486	487	488	489	491	493	494	495	496	497	501	502	503	505	509	520	523	531	535	538	544	546	552	554	555	557	559	562	563	564	567	569	571	572	573	575	576	577	580	581	584	588	591	595	597	598	604	606	607	609	612	613	614	620	621	624	627	629	635	641	642	644	647	650	651	652	655	657	660	662	664	665	669	670	674	676	680	683	689	696	697	699	700	702	704	705	716	718	721	725	727	728	734	735	737	740	741	742	743	747	748	751	753	754	757	764	765	768	774	776	777	783	788	789	790	792	793	796	799	801	805	806	808	809	811	816	819	820	823	824	826	829	830	831	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [722	391	417	408	192	72	57	739	194	715	590	838	775	116	325	808	531	736	396	12	510	189	431	591	495	479	352	199	348	9	492	774	77	29	713	811	8	653	113	704	806	281	473	466	822	514	344	111	595	191	513	675	67	764	25	655	654	39	426	735	251	627	467	612	619	781	219	504	378	749	575	436	230	196	678	341	622	75	782	21	682	758	397	91];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq51'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq52'
        newLIndices = [3	6	9	10	11	16	18	20	25	26	28	29	32	33	34	36	38	39	40	46	49	50	51	54	56	57	60	62	65	67	68	69	71	72	74	80	81	84	92	94	97	100	101	102	110	112	113	114	115	116	117	119	121	123	124	126	128	129	137	139	142	149	151	152	156	159	160	164	173	174	175	177	181	185	187	188	189	197	198	202	204	205	207	208	209	210	211	214	215	218	219	221	223	224	226	230	232	237	240	241	243	244	255	258	268	270	272	273	275	276	279	285	287	289	299	306	307	309	310	312	323	325	326	327	329	331	333	334	335	336	339	340	341	346	347	354	358	359	362	366	368	370	373	376	378	379	380	382	387	388	389	390	393	396	398	399	401	405	412	413	415	416	418	419	421	426	427	429	433	441	442	444	445	447	448	454	456	457	465	466	471	476	479	480	484	487	488	490	493	494	496	500	502	506	510	511	512	514	515	516	518	519	520	521	524	528	529	530	534	535	536	537	546	554	556	559	563	568	569	570	571	574	576	580	586	587	588	591	592	593	596	597	601	602	603	608	610	614	616	617	618	619	624	630	635	637	638	641	643	644	645	649	650	652	653	655	657	659	660	664	665	668	671	672	674	675	677	678	680	685	687	688	689	693	694	695	703	704	705	712	716	718	720	723	725	726	735	737	746	747	749	754	756	757	759	761	762	763	764	765	766	767	768	770	771	775	777	779	781	783	784	785	787	789	792	793	795	798	799	802	803	804	805	806	809	812	813	814	815	816	820	827	832	835	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	5	7	8	12	13	14	15	17	19	21	22	23	24	27	30	31	35	37	41	42	43	44	45	47	48	52	53	55	58	59	61	63	64	66	70	73	75	76	77	78	79	82	83	85	86	87	88	89	90	91	93	95	96	98	99	103	104	105	106	107	108	109	111	118	120	122	125	127	130	131	132	133	134	135	136	138	140	141	143	144	145	146	147	148	150	153	154	155	157	158	161	162	163	165	166	167	168	169	170	171	172	176	178	179	180	182	183	184	186	190	191	192	193	194	195	196	199	200	201	203	206	212	213	216	217	220	222	225	227	228	229	231	233	234	235	236	238	239	242	245	246	247	248	249	250	251	252	253	254	256	257	259	260	261	262	263	264	265	266	267	269	271	274	277	278	280	281	282	283	284	286	288	290	291	292	293	294	295	296	297	298	300	301	302	303	304	305	308	311	313	314	315	316	317	318	319	320	321	322	324	328	330	332	337	338	342	343	344	345	348	349	350	351	352	353	355	356	357	360	361	363	364	365	367	369	371	372	374	375	377	381	383	384	385	386	391	392	394	395	397	400	402	403	404	406	407	408	409	410	411	414	417	420	422	423	424	425	428	430	431	432	434	435	436	437	438	439	440	443	446	449	450	451	452	453	455	458	459	460	461	462	463	464	467	468	469	470	472	473	474	475	477	478	481	482	483	485	486	489	491	492	495	497	498	499	501	503	504	505	507	508	509	513	517	522	523	525	526	527	531	532	533	538	539	540	541	542	543	544	545	547	548	549	550	551	552	553	555	557	558	560	561	562	564	565	566	567	572	573	575	577	578	579	581	582	583	584	585	589	590	594	595	598	599	600	604	605	606	607	609	611	612	613	615	620	621	622	623	625	626	627	628	629	631	632	633	634	636	639	640	642	646	647	648	651	654	656	658	661	662	663	666	667	669	670	673	676	679	681	682	683	684	686	690	691	692	696	697	698	699	700	701	702	706	707	708	709	710	711	713	714	715	717	719	721	722	724	727	728	729	730	731	732	733	734	736	738	739	740	741	742	743	744	745	748	750	751	752	753	755	758	760	769	772	773	774	776	778	780	782	786	788	790	791	794	796	797	800	801	807	808	810	811	817	818	819	821	822	823	824	825	826	828	829	830	831	833	834	836	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [802	373	72	48	525	664	576	287	786	432	790	61	172	640	754	645	243	125	696	643	222	187	262	677	670	464	465	232	567	823	357	361	376	225	544	726	729	599	209	552	106	99	152	116	812	59	652	572	732	389	517	701	424	818	766	29	256	760	285	241	95	712	798	258	695	387	476	450	538	23	406	25	634	261	647	188	444	714	37	41	16	516	453	87];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq53'
        newLIndices = [5	9	12	13	14	21	23	26	33	35	37	38	42	44	46	49	52	53	54	62	66	68	69	73	76	77	82	85	88	90	91	92	94	95	98	107	109	113	124	126	130	133	135	136	148	151	152	154	156	157	159	161	164	167	168	171	173	174	185	188	192	201	204	205	211	215	216	221	234	236	237	240	245	250	252	253	254	265	267	273	276	277	280	281	282	283	285	289	291	295	296	298	301	302	304	309	312	318	322	323	325	326	342	347	361	363	365	367	370	371	375	384	387	389	403	412	414	416	418	421	436	438	439	440	442	445	448	449	450	451	454	455	456	462	464	474	479	481	486	491	493	496	500	504	507	508	509	511	518	519	520	521	525	530	532	533	536	541	550	551	554	556	559	560	563	569	570	573	578	589	590	593	595	598	599	608	610	611	621	622	628	634	638	640	645	649	650	652	656	657	660	665	668	673	679	680	681	684	685	687	690	691	692	693	696	701	702	704	709	710	711	713	726	736	738	742	747	753	754	755	756	759	761	765	772	773	774	778	779	780	783	785	790	791	792	798	800	804	806	807	808	810	817	823	827	829	830	832	834	835	836];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	3	4	6	7	8	10	11	15	16	17	18	19	20	22	24	25	27	28	29	30	31	32	34	36	39	40	41	43	45	47	48	50	51	55	56	57	58	59	60	61	63	64	65	67	70	71	72	74	75	78	79	80	81	83	84	86	87	89	93	96	97	99	100	101	102	103	104	105	106	108	110	111	112	114	115	116	117	118	119	120	121	122	123	125	127	128	129	131	132	134	137	138	139	140	141	142	143	144	145	146	147	149	150	153	155	158	160	162	163	165	166	169	170	172	175	176	177	178	179	180	181	182	183	184	186	187	189	190	191	193	194	195	196	197	198	199	200	202	203	206	207	208	209	210	212	213	214	217	218	219	220	222	223	224	225	226	227	228	229	230	231	232	233	235	238	239	241	242	243	244	246	247	248	249	251	255	256	257	258	259	260	261	262	263	264	266	268	269	270	271	272	274	275	278	279	284	286	287	288	290	292	293	294	297	299	300	303	305	306	307	308	310	311	313	314	315	316	317	319	320	321	324	327	328	329	330	331	332	333	334	335	336	337	338	339	340	341	343	344	345	346	348	349	350	351	352	353	354	355	356	357	358	359	360	362	364	366	368	369	372	373	374	376	377	378	379	380	381	382	383	385	386	388	390	391	392	393	394	395	396	397	398	399	400	401	402	404	405	406	407	408	409	410	411	413	415	417	419	420	422	423	424	425	426	427	428	429	430	431	432	433	434	435	437	441	443	444	446	447	452	453	457	458	459	460	461	463	465	466	467	468	469	470	471	472	473	475	476	477	478	480	482	483	484	485	487	488	489	490	492	494	495	497	498	499	501	502	503	505	506	510	512	513	514	515	516	517	522	523	524	526	527	528	529	531	534	535	537	538	539	540	542	543	544	545	546	547	548	549	552	553	555	557	558	561	562	564	565	566	567	568	571	572	574	575	576	577	579	580	581	582	583	584	585	586	587	588	591	592	594	596	597	600	601	602	603	604	605	606	607	609	612	613	614	615	616	617	618	619	620	623	624	625	626	627	629	630	631	632	633	635	636	637	639	641	642	643	644	646	647	648	651	653	654	655	658	659	661	662	663	664	666	667	669	670	671	672	674	675	676	677	678	682	683	686	688	689	694	695	697	698	699	700	703	705	706	707	708	712	714	715	716	717	718	719	720	721	722	723	724	725	727	728	729	730	731	732	733	734	735	737	739	740	741	743	744	745	746	748	749	750	751	752	757	758	760	762	763	764	766	767	768	769	770	771	775	776	777	781	782	784	786	787	788	789	793	794	795	796	797	799	801	802	803	805	809	811	812	813	814	815	816	818	819	820	821	822	824	825	826	828	831	833	837	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [696	244	337	720	513	826	170	688	561	207	394	331	496	661	87	677	692	292	353	469	574	607	619	318	350	778	466	690	224	504	476	813	70	403	420	73	726	709	352	625	119	494	208	355	671	156	241	383	267	631	779	126	187	552	295	763	761	503	794	314	492	766	435	724	558	375	797	685	154	305	762	309	505	811	782	499	83	28	470	431	730	565	501	396];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq54'
        newLIndices = [1	7	8	9	20	28	40	42	43	44	45	60	64	65	68	72	76	88	97	113	115	124	139	141	144	146	151	153	154	161	163	170	172	173	182	186	190	195	196	199	201	203	213	214	215	220	222	225	227	251	255	261	263	271	273	274	277	282	284	287	290	294	301	308	313	321	327	333	334	343	347	358	360	361	362	365	368	388	393	394	395	407	410	415	420	429	430	434	438	443	444	446	452	461	466	486	492	495	496	502	505	506	513	519	521	522	537	539	545	547	549	555	566	569	571	577	578	579	582	584	586	590	596	599	606	610	611	620	639	642	646	655	666	686	688	694	696	699	702	709	719	727	735	738	739	741	743	748	749	750	761	763	764	767	772	773	781	784	786	787	802	806	808	809	815	821	832	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [2	3	4	5	6	10	11	12	13	14	15	16	17	18	19	21	22	23	24	25	26	27	29	30	31	32	33	34	35	36	37	38	39	41	46	47	48	49	50	51	52	53	54	55	56	57	58	59	61	62	63	66	67	69	70	71	73	74	75	77	78	79	80	81	82	83	84	85	86	87	89	90	91	92	93	94	95	96	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	114	116	117	118	119	120	121	122	123	125	126	127	128	129	130	131	132	133	134	135	136	137	138	140	142	143	145	147	148	149	150	152	155	156	157	158	159	160	162	164	165	166	167	168	169	171	174	175	176	177	178	179	180	181	183	184	185	187	188	189	191	192	193	194	197	198	200	202	204	205	206	207	208	209	210	211	212	216	217	218	219	221	223	224	226	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	252	253	254	256	257	258	259	260	262	264	265	266	267	268	269	270	272	275	276	278	279	280	281	283	285	286	288	289	291	292	293	295	296	297	298	299	300	302	303	304	305	306	307	309	310	311	312	314	315	316	317	318	319	320	322	323	324	325	326	328	329	330	331	332	335	336	337	338	339	340	341	342	344	345	346	348	349	350	351	352	353	354	355	356	357	359	363	364	366	367	369	370	371	372	373	374	375	376	377	378	379	380	381	382	383	384	385	386	387	389	390	391	392	396	397	398	399	400	401	402	403	404	405	406	408	409	411	412	413	414	416	417	418	419	421	422	423	424	425	426	427	428	431	432	433	435	436	437	439	440	441	442	445	447	448	449	450	451	453	454	455	456	457	458	459	460	462	463	464	465	467	468	469	470	471	472	473	474	475	476	477	478	479	480	481	482	483	484	485	487	488	489	490	491	493	494	497	498	499	500	501	503	504	507	508	509	510	511	512	514	515	516	517	518	520	523	524	525	526	527	528	529	530	531	532	533	534	535	536	538	540	541	542	543	544	546	548	550	551	552	553	554	556	557	558	559	560	561	562	563	564	565	567	568	570	572	573	574	575	576	580	581	583	585	587	588	589	591	592	593	594	595	597	598	600	601	602	603	604	605	607	608	609	612	613	614	615	616	617	618	619	621	622	623	624	625	626	627	628	629	630	631	632	633	634	635	636	637	638	640	641	643	644	645	647	648	649	650	651	652	653	654	656	657	658	659	660	661	662	663	664	665	667	668	669	670	671	672	673	674	675	676	677	678	679	680	681	682	683	684	685	687	689	690	691	692	693	695	697	698	700	701	703	704	705	706	707	708	710	711	712	713	714	715	716	717	718	720	721	722	723	724	725	726	728	729	730	731	732	733	734	736	737	740	742	744	745	746	747	751	752	753	754	755	756	757	758	759	760	762	765	766	768	769	770	771	774	775	776	777	778	779	780	782	783	785	788	789	790	791	792	793	794	795	796	797	798	799	800	801	803	804	805	807	810	811	812	813	814	816	817	818	819	820	822	823	824	825	826	827	828	829	830	831	833	834	835	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [239	693	160	370	329	689	564	173	265	111	556	473	141	122	393	748	454	28	45	660	835	313	645	297	434	579	708	267	527	789	62	474	334	249	213	610	799	150	625	157	792	640	338	580	396	642	283	58	467	719	153	341	589	31	742	598	437	144	389	404	774	665	747	810	817	723	371	179	306	543	429	581	763	737	409	736	89	40	232	441	403	683	764	327];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq55'
        newLIndices = [1	7	8	20	28	40	44	72	97	115	124	141	146	151	153	154	182	190	195	196	199	215	222	227	251	261	263	271	273	274	287	301	308	327	333	334	343	347	362	394	395	410	430	434	443	444	446	452	461	486	492	495	502	505	513	521	522	537	539	545	549	555	569	582	584	586	590	599	639	642	655	696	699	709	727	738	739	741	748	749	750	767	786	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [2	3	4	5	6	9	10	11	12	13	14	15	16	17	18	19	21	22	23	24	25	26	27	29	30	31	32	33	34	35	36	37	38	39	41	42	43	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	116	117	118	119	120	121	122	123	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	142	143	144	145	147	148	149	150	152	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	183	184	185	186	187	188	189	191	192	193	194	197	198	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	216	217	218	219	220	221	223	224	225	226	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	252	253	254	255	256	257	258	259	260	262	264	265	266	267	268	269	270	272	275	276	277	278	279	280	281	282	283	284	285	286	288	289	290	291	292	293	294	295	296	297	298	299	300	302	303	304	305	306	307	309	310	311	312	313	314	315	316	317	318	319	320	321	322	323	324	325	326	328	329	330	331	332	335	336	337	338	339	340	341	342	344	345	346	348	349	350	351	352	353	354	355	356	357	358	359	360	361	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377	378	379	380	381	382	383	384	385	386	387	388	389	390	391	392	393	396	397	398	399	400	401	402	403	404	405	406	407	408	409	411	412	413	414	415	416	417	418	419	420	421	422	423	424	425	426	427	428	429	431	432	433	435	436	437	438	439	440	441	442	445	447	448	449	450	451	453	454	455	456	457	458	459	460	462	463	464	465	466	467	468	469	470	471	472	473	474	475	476	477	478	479	480	481	482	483	484	485	487	488	489	490	491	493	494	496	497	498	499	500	501	503	504	506	507	508	509	510	511	512	514	515	516	517	518	519	520	523	524	525	526	527	528	529	530	531	532	533	534	535	536	538	540	541	542	543	544	546	547	548	550	551	552	553	554	556	557	558	559	560	561	562	563	564	565	566	567	568	570	571	572	573	574	575	576	577	578	579	580	581	583	585	587	588	589	591	592	593	594	595	596	597	598	600	601	602	603	604	605	606	607	608	609	610	611	612	613	614	615	616	617	618	619	620	621	622	623	624	625	626	627	628	629	630	631	632	633	634	635	636	637	638	640	641	643	644	645	646	647	648	649	650	651	652	653	654	656	657	658	659	660	661	662	663	664	665	666	667	668	669	670	671	672	673	674	675	676	677	678	679	680	681	682	683	684	685	686	687	688	689	690	691	692	693	694	695	697	698	700	701	702	703	704	705	706	707	708	710	711	712	713	714	715	716	717	718	719	720	721	722	723	724	725	726	728	729	730	731	732	733	734	735	736	737	740	742	743	744	745	746	747	751	752	753	754	755	756	757	758	759	760	761	762	763	764	765	766	768	769	770	771	772	773	774	775	776	777	778	779	780	781	782	783	784	785	787	788	789	790	791	792	793	794	795	796	797	798	799	800	801	802	803	804	805	806	807	808	809	810	811	812	813	814	815	816	817	818	819	820	821	822	823	824	825	826	827	828	829	830	831	832	833	834	835	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [314	487	98	49	818	238	496	800	155	161	283	772	323	226	126	327	309	108	357	75	503	9	468	644	192	365	463	50	402	520	179	676	783	682	407	224	599	190	766	832	480	138	72	203	805	723	555	574	182	455	639	318	777	71	252	401	48	567	435	413	646	668	613	247	351	582	85	827	208	404	747	545	239	223	650	696	488	195	68	637	444	718	47	442];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq56'
        newMIndices = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	259	260	261	262	263	264	265	266	267	268	269	270	271	272	273	274	275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303	304	305	306	307	308	309	310	311	312	313	314	315	316	317	318	319	320	321	322	323	324	325	326	327	328	329	330	331	332	333	334	335	336	337	338	339	340	341	342	343	344	345	346	347	348	349	350	351	352	353	354	355	356	357	358	359	360	361	362	363	364	365	366	367	368	369	370	371	372	373	374	375	376	377	378	379	380	381	382	383	384	385	386	387	388	389	390	391	392	393	394	395	396	397	398	399	400	401	402	403	404	405	406	407	408	409	410	411	412	413	414	415	416	417	418	419	420	421	422	423	424	425	426	427	428	429	430	431	432	433	434	435	436	437	438	439	440	441	442	443	444	445	446	447	448	449	450	451	452	453	454	455	456	457	458	459	460	461	462	463	464	465	466	467	468	469	470	471	472	473	474	475	476	477	478	479	480	481	482	483	484	485	486	487	488	489	490	491	492	493	494	495	496	497	498	499	500	501	502	503	504	505	506	507	508	509	510	511	512	513	514	515	516	517	518	519	520	521	522	523	524	525	526	527	528	529	530	531	532	533	534	535	536	537	538	539	540	541	542	543	544	545	546	547	548	549	550	551	552	553	554	555	556	557	558	559	560	561	562	563	564	565	566	567	568	569	570	571	572	573	574	575	576	577	578	579	580	581	582	583	584	585	586	587	588	589	590	591	592	593	594	595	596	597	598	599	600	601	602	603	604	605	606	607	608	609	610	611	612	613	614	615	616	617	618	619	620	621	622	623	624	625	626	627	628	629	630	631	632	633	634	635	636	637	638	639	640	641	642	643	644	645	646	647	648	649	650	651	652	653	654	655	656	657	658	659	660	661	662	663	664	665	666	667	668	669	670	671	672	673	674	675	676	677	678	679	680	681	682	683	684	685	686	687	688	689	690	691	692	693	694	695	696	697	698	699	700	701	702	703	704	705	706	707	708	709	710	711	712	713	714	715	716	717	718	719	720	721	722	723	724	725	726	727	728	729	730	731	732	733	734	735	736	737	738	739	740	741	742	743	744	745	746	747	748	749	750	751	752	753	754	755	756	757	758	759	760	761	762	763	764	765	766	767	768	769	770	771	772	773	774	775	776	777	778	779	780	781	782	783	784	785	786	787	788	789	790	791	792	793	794	795	796	797	798	799	800	801	802	803	804	805	806	807	808	809	810	811	812	813	814	815	816	817	818	819	820	821	822	823	824	825	826	827	828	829	830	831	832	833	834	835	836	837	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [639	7	569	590	539	461	182	642	190	308	738	709	333	263	502	750	748	486	273	699	362	739	28	434	584	146	274	153	261	327	444	40	446	222	195	196	124	767	749	655	582	141	287	151	1	251	555	495	430	347	227	395	599	786	452	586	505	97	394	271	72	115	154	522	334	537	199	8	410	215	727	696	301	20	513	838	741	44	343	443	521	545	492	549];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()




        % Quad Seq 8 - 30 arcmin
    case 'quadSeq57'

        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()



        % Quad Seq 10 - 30 arcmin
    case 'quadSeq58' % Blue print
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq59'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq60'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq61'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq62'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq63'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq64'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq65'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq66'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newLIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()

    case 'quadSeq67'
        newLIndices = [3	5	7	8	9	13	14	16	20	21	22	23	25	26	27	29	31	32	33	37	39	40	41	43	45	46	49	51	53	55	56	57	58	59	61	66	67	69	75	76	78	80	81	82	88	90	91	92	93	94	95	96	98	100	101	103	105	106	112	114	117	122	124	125	128	130	131	134	141	142	143	145	148	151	152	153	154	160	161	164	166	167	169	170	171	172	173	175	176	178	179	180	182	183	184	187	189	192	194	195	196	197	205	208	215	216	217	218	220	221	223	228	230	231	239	244	245	246	247	249	257	259	260	261	263	265	267	268	269	270	272	273	274	277	278	283	286	287	290	293	294	296	298	300	302	303	304	305	309	310	311	312	314	317	318	319	321	324	329	330	332	333	335	336	338	342	343	345	348	354	355	357	358	360	361	366	367	368	374	375	379	383	385	386	389	391	392	393	395	396	398	401	403	406	409	410	411	413	414	415	417	418	419	420	422	425	426	427	430	431	432	433	440	446	448	450	453	457	458	459	460	462	464	467	472	473	474	477	478	479	481	482	485	486	487	491	493	496	497	498	499	500	504	509	513	515	516	518	520	521	522	525	526	528	529	531	532	534	535	538	539	541	543	544	546	547	549	550	551	555	556	557	558	561	562	563	569	570	571	576	579	581	583	585	587	588	595	597	604	605	607	612	614	615	616	618	619	620	621	622	623	624	625	627	628	631	632	634	636	638	639	640	641	643	645	646	647	650	651	653	654	655	656	657	659	661	662	663	664	665	667	670	673	675	677	682	683	685	687	691	695	696	700	701	704	705	706	709	711	713	714	715	717	720	721	722	723	725	728	732	734	735	740	741	742	743	744	747	749	751	752	755	758	760	761	762	764	766	767	768	769	770	773	775	776	778	779	784	785	787	789	790	792	794	798	799	800	803	804	805	807	809	810	811	814	815	817	818	819	820	821	823	824	825	827	828	829	835	838];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = [1	2	4	6	10	11	12	15	17	18	19	24	28	30	34	35	36	38	42	44	47	48	50	52	54	60	62	63	64	65	68	70	71	72	73	74	77	79	83	84	85	86	87	89	97	99	102	104	107	108	109	110	111	113	115	116	118	119	120	121	123	126	127	129	132	133	135	136	137	138	139	140	144	146	147	149	150	155	156	157	158	159	162	163	165	168	174	177	181	185	186	188	190	191	193	198	199	200	201	202	203	204	206	207	209	210	211	212	213	214	219	222	224	225	226	227	229	232	233	234	235	236	237	238	240	241	242	243	248	250	251	252	253	254	255	256	258	262	264	266	271	275	276	279	280	281	282	284	285	288	289	291	292	295	297	299	301	306	307	308	313	315	316	320	322	323	325	326	327	328	331	334	337	339	340	341	344	346	347	349	350	351	352	353	356	359	362	363	364	365	369	370	371	372	373	376	377	378	380	381	382	384	387	388	390	394	397	399	400	402	404	405	407	408	412	416	421	423	424	428	429	434	435	436	437	438	439	441	442	443	444	445	447	449	451	452	454	455	456	461	463	465	466	468	469	470	471	475	476	480	483	484	488	489	490	492	494	495	501	502	503	505	506	507	508	510	511	512	514	517	519	523	524	527	530	533	536	537	540	542	545	548	552	553	554	559	560	564	565	566	567	568	572	573	574	575	577	578	580	582	584	586	589	590	591	592	593	594	596	598	599	600	601	602	603	606	608	609	610	611	613	617	626	629	630	633	635	637	642	644	648	649	652	658	660	666	668	669	671	672	674	676	678	679	680	681	684	686	688	689	690	692	693	694	697	698	699	702	703	707	708	710	712	716	718	719	724	726	727	729	730	731	733	736	737	738	739	745	746	748	750	753	754	756	757	759	763	765	771	772	774	777	780	781	782	783	786	788	791	793	795	796	797	801	802	806	808	812	813	816	822	826	830	831	832	833	834	836	837];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [699	670	768	115	421	338	145	478	504	178	431	819	405	574	340	29	241	658	285	69	419	300	604	428	655	665	154	101	666	517	14	724	416	439	488	611	686	307	68	587	265	837	296	659	141	103	697	35	543	579	345	299	771	314	346	123	255	246	698	193	242	318	550	112	675	65	357	24	580	539	165	522	427	651	766	689	320	273	372	195	705	354	779	326];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % CenterPatch
        newMIndices = [390 414 379 431 371 415 403 362 445 404 361 370 430 ...
            350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);
        theConeMosaic.visualizeMosaic()










    case 'quadSeq68'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [390 414 415 403 362 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430	379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq69'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [414 415 403 362 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430 379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq70'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [415 403 362 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430 379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq71'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [403 362 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 415];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430	379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq72'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [362 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 415 403];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430	379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq73'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 415 403 362];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430	379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq74'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 415 403 362 404];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430	379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq75'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 415 403 362 404 361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430	379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

        theConeMosaic.visualizeMosaic(); hold on;

    case 'quadSeq76'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 415 403 362 404 361 370];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [430 379];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);



        % 77-87 Repeat of 57-67 with surround 30% L
    case 'quadSeq77'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; %
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq78'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq79'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq80'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};  % [403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq81'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq82'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq83'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq84'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq85'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq86'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq87'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; %
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);





        % 88 - 98 Repeat of 57-67 with surround 70%L
    case 'quadSeq88'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; %
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq89'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq90'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq91'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq92'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq93'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq94'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq95'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq96'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq97'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq98'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; %
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);







        % 99 - 109 Repeat of 57-67 with surround alt L/S (QS34)
    case 'quadSeq99'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);
        %
        %         newMIndices = storedQuadIndices{1,seqNum}{2}; %
        %         theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq100'
        seqNum = str2double(chrom(8:end));

        newLIndices =  storedQuadIndices{1,seqNum}{1}; % [379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices =  storedQuadIndices{1,seqNum}{2}; % [390 414]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq101'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq102'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq103'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq104'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq105'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq106'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq107'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq108'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; % [416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469]
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq109'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1}; %
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2}; % [390 414 379 431 371 415 403 362 445 404 361 370 430 350 459 429 349 469 416 369 351];
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3}; % [404	361];
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);






        % OVERHAUL, RESET SO AT NEW 10 ARCMIN PROPORTIONS

    case 'quadSeq110'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq111'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq112'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq113'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq114'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq115'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq116'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq117'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq118'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);







        % 3.5 with proper prop no fixed S
    case 'quadSeq119'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq120'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq121'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq122'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq123'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq124'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq125'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq126'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq127'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);







        % 3.5 proper prop with fixes S
    case 'quadSeq128'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq129'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq130'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq131'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq132' % The 50% case
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq133'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);


    case 'quadSeq134'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq135'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq136'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);









    case 'quadSeq137'    % 3.5 position sanity checks
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq138'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq139'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq140'
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq141'    % 3.5 all L annulus surround
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);

    case 'quadSeq142'   % 3.5 all M annulus surround
        seqNum = str2double(chrom(8:end));

        newLIndices = storedQuadIndices{1,seqNum}{1};
        theConeMosaic.Mosaic.reassignTypeOfCones(newLIndices, cMosaic.LCONE_ID);

        newMIndices = storedQuadIndices{1,seqNum}{2};
        theConeMosaic.Mosaic.reassignTypeOfCones(newMIndices, cMosaic.MCONE_ID);

        newSIndices = storedQuadIndices{1,seqNum}{3};
        theConeMosaic.Mosaic.reassignTypeOfCones(newSIndices, cMosaic.SCONE_ID);














    case 'quadSeqNew'
        open tempConeProps.m
        keyboard()
    otherwise
        error(['Entered chrom has no alloted override sequence']);
end

% Save the updated cell of QuadSeqs for future
save(fullfile(storedDir, "storedQuadIndices.mat"), "storedQuadIndices")

end


%% Useful tools for building QuadSeqs

% % Visaulize index numbers over each cone in the mosaic
%
% theConeMosaic.visualizeMosaic()
% for i=1:length(theConeMosaic.Mosaic.coneTypes)
%     txt = int2str(i);
%     t = text((theConeMosaic.Mosaic.coneRFpositionsDegs(i,1)-0.005), theConeMosaic.Mosaic.coneRFpositionsDegs(i,2),txt);
%     t.FontSize=11;
%     hold on
% end

% % Call function to manually select and swap cones (see mouseClick)
%
% g = @(x, y) mouseClick(x, y, theConeMosaic);
% set(gcf, 'WindowButtonDownFcn', g)
% keyboard

% Save the new sequence selections to the stored file. If using mouseClick,
% can skip to the last line
%
% bookKeep(1) = {theConeMosaic.Mosaic.lConeIndices'};
% bookKeep(2) = {theConeMosaic.Mosaic.mConeIndices'};
% bookKeep(3) = {theConeMosaic.Mosaic.sConeIndices'};
% storedQuadIndices(1,(seqNum)) = {bookKeep};


