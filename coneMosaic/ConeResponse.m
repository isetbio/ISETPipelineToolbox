classdef ConeResponse < handle
    %CONERESPONSE Class that computes cone responses from RGB image.
    %
    % Syntax: retina = ConeResponse([varargin]);
    %
    % Description:
    %   This wrapper class can compute the mean cone excitations and the
    %   optical image for a given retinal mosaic and RGB image, use
    %   function provided by isetbio.
    %
    %   Although these could be varied in the future with key/value pairs,
    %   this object uses:
    %     'whichOptics'       - 'wvfHuman'  Wavefront based human eye and optics.
    %                            Inherits isebio default wavefront optics and
    %                            pupil size.
    %     'whichMosaic'       - 'coneMosaicHex' Mosaic based on isetbio
    %                            coneMosaicHex method, as specified in detail
    %                            by additional key value pairs.
    %
    % ConeResponse Properties:
    %   Display        - Display used for the image (isetbio/displayCreate)
    %   FovealDegree   - Degree of visual field we are modeling
    %   PSF            - Point spread function of the eye (isetbio/oiCreate)
    %   Mosaic         - Photoreceptor mosaic (isetbio/coneMosaicHex)
    %   LastResponse   - Recent mosaic excitation computed
    %   LastOI         - Recent optical image computed
    %
    % ConeResponse Methods:
    %   ConeResponse         - Constructor for ConeResponse Object
    %   compute              - Compute mosaic excitation given a RGB image
    %   visualizeMosaic      - Visualize the cone mosaic being simulated
    %   visualizeExcitation  - Visualize the mosaic excitation pattern
    %   visualizeOI          - Visualize the optical image
    %
    % Inputs:
    %   None.
    %
    % Outputs:
    %   ConeResponse Object.
    %
    % Optional key/value pairs:
    %   'fovealDegree'        - Double. The degree of visual angle we are using
    %                           for our cone mosaic.
    %                           Default is 1.0 (degree).
    %   'eccBasedConeDensity' - Boolean. Vary cone density based on
    %                           eccentricity or not.
    %                           Default is False.
    %   'eccBasedConeQuantal' - Boolean. Vary cone quantal efficiency based
    %                           on eccentricity or not.
    %                           Default is False.
    %   'display'             - Display object used in isetbio routine.
    %                           Default is 'LCD-Apple'.
    %   'viewDistance'        - Double. View distance from the screen.
    %                           Default is 0.57.
    %   'spatialDensity'      - Vector. Density of L, M and S cone.
    %                           Default is [].
    
    
    properties (Access = public)
        Display;
        FovealDegree;
        PSF;
        diffLmt;
        Mosaic;
        LastResponse;
        LastOI;
        PupilSize;
        DefaultMosaic;
        
        L_Cone_Idx = 2;
        M_Cone_Idx = 3;
        S_Cone_Idx = 4;
    end
    
    methods(Static)
        function oiDiffLmt = psfDiffLmt(pupilDiamMM,varargin)

            % Parse input
            p = inputParser;
            p.addParameter('defocusDiopters', 0, @isnumeric);
            parse(p, varargin{:});

            if ~exist('pupilSize', 'var')
                pupilDiamMM = 3.0;
            end
            
            % Set up wavefront optics object
            wave = (400:10:700)';
            accommodatedWavelength = 530;
            pupilDiameterMm = pupilDiamMM;
            
            % Zero Zernike coefficients (diffraction limited)
            zCoeffs = zeros(66,1);
            wvfP = wvfCreate('calc wavelengths', wave, 'zcoeffs', zCoeffs, ...
                'name', sprintf('human-%d', pupilDiameterMm));
            wvfP = wvfSet(wvfP, 'measured pupil size', pupilDiameterMm);
            wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMm);
            wvfP = wvfSet(wvfP, 'measured wavelength', accommodatedWavelength);

            % Add in specified defocus
            defocusMicrons = wvfDefocusDioptersToMicrons(p.Results.defocusDiopters,pupilDiamMM);
            wvfP = wvfSet(wvfP,'zcoeffs', defocusMicrons, 'defocus');
            
            % Compute pupil function using 'no lca' key/value pair to turn off LCA.
            wvfPNoLca = wvfComputePupilFunction(wvfP,false,'no lca',true);
            wvfPNoLca = wvfComputePSF(wvfPNoLca);
            
            % Optics object
            oiDiffLmt = wvf2oi(wvfPNoLca);
            opticsNoLca = oiGet(oiDiffLmt, 'optics');
            opticsNoLca = opticsSet(opticsNoLca, 'model', 'shift invariant');
            opticsNoLca = opticsSet(opticsNoLca, 'name', 'human-wvf-nolca');
            oiDiffLmt = oiSet(oiDiffLmt,'optics',opticsNoLca);
        end
        
        function psfNoLCA = psfNoLCA(pupilSize)
            if ~exist('pupilSize', 'var')
                pupilSize = 3.0;
            end
            
            % Set up wavefront optics object
            wave = (400:10:700)';
            accommodatedWavelength = 530;
            pupilDiameterMm = pupilSize;
            
            % get Zernike coefficients for humans
            zCoeffs = wvfLoadThibosVirtualEyes(pupilDiameterMm);
            wvfP = wvfCreate('calc wavelengths', wave, 'zcoeffs', zCoeffs, ...
                'name', sprintf('human-%d', pupilDiameterMm));
            wvfP = wvfSet(wvfP, 'measured pupil size', pupilDiameterMm);
            wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMm);
            wvfP = wvfSet(wvfP, 'measured wavelength', accommodatedWavelength);
            
            % Compute pupil function using 'no lca' key/value pair to turn off LCA.
            wvfPNoLca = wvfComputePupilFunction(wvfP,false,'no lca',true);
            wvfPNoLca = wvfComputePSF(wvfPNoLca);
            
            % Optics object
            psfDiffLmt = wvf2oi(wvfPNoLca);
            opticsNoLca = oiGet(psfDiffLmt, 'optics');
            opticsNoLca = opticsSet(opticsNoLca, 'model', 'shift invariant');
            opticsNoLca = opticsSet(opticsNoLca, 'name', 'human-wvf-nolca');
            psfNoLCA = oiSet(psfDiffLmt,'optics',opticsNoLca);
        end
        
    end
    
    methods (Access = public)
        
        function this = ConeResponse(varargin)
            % CONERESPONSE  Construct ConeResponse object.
            %   Construct ConeResponse object with optional argument:
            %   FovealDegree, eccBasedConeDensity, eccBasedConeQuantal.
            display = displayCreate('CRT12BitDisplay');
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('fovealDegree', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('eccBasedConeDensity', false, @islogical);
            p.addParameter('eccBasedConeQuantal', false, @islogical);
            p.addParameter('display', display);
            p.addParameter('viewDistance', 0.57, @(x) (isnumeric(x) && numel(x) == 1));
            p.addParameter('spatialDensity', []);
            p.addParameter('pupilSize', 3.0, @(x) (isnumeric(x) && numel(x) == 1));
            p.addParameter('override', false, @islogical);
            p.addParameter('integrationTime', 0.2, @(x) (isnumeric(x) && numel(x) == 1));
            
            parse(p, varargin{:});
            this.FovealDegree   = p.Results.fovealDegree;
            spatialDensity      = p.Results.spatialDensity;
            eccBasedConeDensity = p.Results.eccBasedConeDensity;
            eccBasedConeQuantal = p.Results.eccBasedConeQuantal;
            integrationTime     = p.Results.integrationTime;
            
            if(~p.Results.override)
                % Create cone moasic
                fprintf('Create cone moasic object: \n');
                
                if(isempty(spatialDensity))
                    theMosaic = coneMosaicHex(5, ...                               % hex lattice sampling factor
                        'fovDegs', this.FovealDegree, ...                          % match mosaic width to stimulus size
                        'eccBasedConeDensity', eccBasedConeDensity, ...            % cone density varies with eccentricity
                        'eccBasedConeQuantalEfficiency', eccBasedConeQuantal, ...  % cone quantal efficiency varies with eccentricity
                        'integrationTime', integrationTime, ...                                % 0.2s default integration time
                        'maxGridAdjustmentIterations', 50);                        % terminate iterative lattice adjustment after 50 iterations
                else
                    theMosaic = coneMosaicHex(5, ...
                        'fovDegs', this.FovealDegree, ...
                        'spatialDensity', spatialDensity, ...
                        'sConeMinDistanceFactor', 0, ...
                        'sConeFreeRadiusMicrons', 0, ...
                        'integrationTime', integrationTime, ...
                        'maxGridAdjustmentIterations', 50);
                end
                
                % Poisson noise model, mean response
                theMosaic.noiseFlag = 'none';
                this.Mosaic = theMosaic;
                this.DefaultMosaic = this.Mosaic.pattern;
                
                % Point spread function of human eye
                this.PupilSize = p.Results.pupilSize;
                this.PSF = oiCreate('human', p.Results.pupilSize);
            else
                fprintf('Override, will not create mosaic and optics object \n');
            end
            
            % Display
            this.Display = p.Results.display;
            this.Display.dist = p.Results.viewDistance;
            this.diffLmt = this.psfDiffLmt();
        end
        
        function visualizeMosaic(this)
            this.Mosaic.visualizeGrid(...
                'backgroundColor', [1 1 1], ...
                'ticksInVisualDegs', true, ...
                'generateNewFigure', true);
            set(gca, 'XTick', [], 'YTick', []);
            set(gca,'YDir','reverse');
        end
        
        function [excitation, theOI, linearizedImage, allCone, L, M, S] = compute(this, image)
            % COMPUTE   Compute optical image and cone mosaic excitation.
            %
            % Syntax:
            %   [excitation, OI, allCone, L, M, S] = this.compute(image)
            %
            % Description:
            %   Compute OI and cone excitation array of provided RGB image
            %   input, can also return L, M, S cone responses separately.
            %
            % Inputs
            %   image        - Input RGB image
            %
            % Outputs:
            %   excitation   - N by N array representing the mean excitation
            %                  of cone mosaic
            %   theOI        - Struct of the optical image on the retina
            %   allCone      - Vector representing the excitation of all cones
            %   L            - Vector of L cone excitation
            %   M            - Vector of M cone excitation
            %   S            - Vector of S cone excitation
            %
            % Optional key/value pairs:
            %   None.
            
            meanLuminanceCdPerM2 = [];
            [realizedStimulusScene, ~, linearizedImage] = sceneFromFile(image, 'rgb', ...
                meanLuminanceCdPerM2, this.Display);
            
            % set the angular scene width
            realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', this.FovealDegree);
            
            % optics
            theOI = oiCompute(this.PSF, realizedStimulusScene);
            this.LastOI = theOI;
            
            % cone excitations
            nTrialsNum = 2;
            emPath = zeros(nTrialsNum, 1, 2);
            
            % compute mosaic excitation responses
            % without the eye movement path
            this.LastResponse = this.Mosaic.compute(theOI, 'emPath', emPath);
            
            sizeExci   = size(this.LastResponse);
            excitation = reshape(this.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            
            % LMS cone excitation
            L = this.getConetypeResponse(this.L_Cone_Idx);
            M = this.getConetypeResponse(this.M_Cone_Idx);
            S = this.getConetypeResponse(this.S_Cone_Idx);
            
            allCone = [L; M; S];
            
        end
        
        function [excitation, opticalImage, linearizedImage, allCone, L, M, S] = computeDiffLmt(this, image)
            meanLuminanceCdPerM2 = [];
            [realizedStimulusScene, ~, linearizedImage] = sceneFromFile(image, 'rgb', ...
                meanLuminanceCdPerM2, this.Display);
            
            % set the angular scene width
            realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', this.FovealDegree);
            
            opticalImage = oiCompute(this.diffLmt, realizedStimulusScene);
            [excitation, allCone, L, M, S] = this.computeWithOI(opticalImage);
        end
        
        function [excitation, allCone, L, M, S] = computeWithOI(this, opticalImage)
            this.LastOI = opticalImage;
            
            % cone excitations
            nTrialsNum = 2;
            emPath = zeros(nTrialsNum, 1, 2);
            
            % compute mosaic excitation responses
            % without the eye movement path
            this.LastResponse = this.Mosaic.compute(opticalImage, 'emPath', emPath);
            
            sizeExci   = size(this.LastResponse);
            excitation = reshape(this.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            
            % LMS cone excitation
            L = this.getConetypeResponse(this.L_Cone_Idx);
            M = this.getConetypeResponse(this.M_Cone_Idx);
            S = this.getConetypeResponse(this.S_Cone_Idx);
            
            allCone = [L; M; S];
        end
        
        function [excitation, allCone, L, M, S] = computeWithScene(this, inputScene)
            opticalImage = oiCompute(this.PSF, inputScene);
            [excitation, allCone, L, M, S] = this.computeWithOI(opticalImage);
        end
        
        function allCone = computeHyperspectral(this, wave, scale, image)
            scene = sceneCreate('whitenoise');
            
            scene.spectrum.wave = wave;
            scene.data.photons = image .* scale;
            scene = sceneSet(scene, 'fov', this.FovealDegree);
            
            opticalImage = oiCompute(scene, this.PSF);
            [~, allCone] = this.computeWithOI(opticalImage);
        end
        
        function [excitation, allCone, L, M, S] = computeWithSceneDiffLmt(this, inputScene)
            opticalImage = oiCompute(this.diffLmt, inputScene);
            [excitation, allCone, L, M, S] = this.computeWithOI(opticalImage);
        end
        
        function [allCone, coneCount] = coneExcitationRnd(this, factor, type)
            nCone = sum(sum(this.Mosaic.pattern ~= 0));
            sizeExci   = size(this.LastResponse);
            excitation = reshape(this.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            nghSize    = ceil(sqrt(sizeExci(2) * sizeExci(3) / nCone));
            
            switch type
                case 'L'
                    coneType = this.L_Cone_Idx;
                case 'M'
                    coneType = this.M_Cone_Idx;
                case 'S'
                    coneType = this.S_Cone_Idx;
                otherwise
                    error('Unknown cone type.');
            end
            coneIdx = find(this.Mosaic.pattern == coneType);
            exciIdx = datasample(coneIdx, 1);
            excitation(exciIdx) = excitation(exciIdx) * factor;
            
            J = floor(exciIdx/sizeExci(3)); I = exciIdx - sizeExci(3) * J; J = J + 1;
            nghMtx = this.Mosaic.pattern(max(1, I - nghSize):min(sizeExci(2), I + nghSize), max(1, J - nghSize):min(sizeExci(3), J + nghSize));
            coneCount = [sum(sum(nghMtx == this.L_Cone_Idx)), sum(sum(nghMtx == this.M_Cone_Idx)), sum(sum(nghMtx == this.S_Cone_Idx))];
            
            this.LastResponse(1, :, :) = excitation;
            this.LastResponse(2, :, :) = excitation;
            
            L = this.getConetypeResponse(this.L_Cone_Idx);
            M = this.getConetypeResponse(this.M_Cone_Idx);
            S = this.getConetypeResponse(this.S_Cone_Idx);
            
            allCone = [L; M; S];
        end
        
        function mosaic = getMosaic(this)
            % Return the cone mosaic object
            mosaic = this.Mosaic;
        end
        
        function visualizeExcitation(this, asSubplot)
            if ~exist('asSubplot','var')
                % third parameter does not exist, so default it to something
                asSubplot = false;
            end
            if ~asSubplot
                figure();
            end
            this.Mosaic.renderActivationMap(gca, squeeze(this.LastResponse(1,:,:)), ...
                'mapType', 'modulated disks', ...
                'showColorBar', ~asSubplot, ...
                'labelColorBarTicks', ~asSubplot, ...
                'titleForColorBar', 'R*/cone/tau', ...
                'showXLabel', ~asSubplot, ...
                'showYLabel', ~asSubplot, ...
                'showXTicks', ~asSubplot, ...
                'showYTicks', ~asSubplot);
            set(gca, 'XTick', [], 'YTick', []);
            set(gca,'YDir','reverse');
            
            if(~asSubplot)
                title('Cone Excitation Pattern');
            end
        end
        
        function visualizeOI(this)
            % Visualize optical image
            visualizeOpticalImage(this.LastOI, 'displayRadianceMaps', false, ...
                'displayRetinalContrastProfiles', false);
            set(gca,'YDir','reverse');
        end
        
        function rgbOI = rgbOpticalImage(this)
            rgbOI = oiGet(this.LastOI, 'rgb image');
        end

        function setConeRatio(this, l, m)
            assert(l + m <= 1.0);
            [~, ~, ~, numCone] = this.coneCount();

            numL = floor(numCone * l);
            numM = floor(numCone * m);
            numS = numCone - numL - numM;

            typeArray = [ones(1, numL) * this.L_Cone_Idx, ... 
                         ones(1, numM) * this.M_Cone_Idx, ...
                         ones(1, numS) * this.S_Cone_Idx];

            typeArray = typeArray(randperm(length(typeArray)));
            coneIdx = find(this.Mosaic.pattern > 0);
            assert(length(typeArray) == length(coneIdx));

            for idx = 1 : length(coneIdx)
                this.Mosaic.pattern(coneIdx(idx)) = typeArray(idx);
            end
        end

        
        function [L, M, S, numCone] = coneCount(this)
            L = sum(this.Mosaic.pattern == this.L_Cone_Idx, 'all');
            M = sum(this.Mosaic.pattern == this.M_Cone_Idx, 'all');
            S = sum(this.Mosaic.pattern == this.S_Cone_Idx, 'all');
            numCone = L + M + S;
        end
        
        function resetCone(this)
            this.Mosaic.pattern = this.DefaultMosaic;
        end
        
        function resetSCone(this)
            coneIdx = find(this.Mosaic.pattern == this.S_Cone_Idx);
            for idx = 1:length(coneIdx)
                if rand < 0.5
                    this.Mosaic.pattern(coneIdx(idx)) = this.L_Cone_Idx;
                else
                    this.Mosaic.pattern(coneIdx(idx)) = this.M_Cone_Idx;
                end
            end
        end
        
        function reassignSCone(this, ratio)
            [~, ~, S, numCone] = this.coneCount();
            deltaS = floor(ratio * numCone) - S;
            
            coneIdx = find(this.Mosaic.pattern == this.L_Cone_Idx | this.Mosaic.pattern == this.M_Cone_Idx);
            reassignIdx = sort(datasample(1:length(coneIdx), deltaS, 2, 'Replace', false));
            reassignIdx = coneIdx(reassignIdx);
            
            this.Mosaic.pattern(reassignIdx) = this.S_Cone_Idx;
        end
        
        function reassignCone(this, ratio, target, replace, showMosaic)
            if ~exist('showMosaic', 'var')
                showMosaic = true;
            end
            
            [L, M, S, numCone] = this.coneCount();
            coneCount = [L, M, S];
            
            numTarget = floor(ratio * numCone);
            if(numTarget > coneCount(target - 1))
                error('Targe cone ratio is higher than the original mosaic.');
            end
            
            numReassign = coneCount(target - 1) - numTarget;
            coneIdx = find(this.Mosaic.pattern == target);
            reassignIdx = sort(datasample(1:length(coneIdx), numReassign, 2, 'Replace', false));
            reassignIdx = coneIdx(reassignIdx);
            
            this.Mosaic.pattern(reassignIdx) = replace;
            if showMosaic
                this.visualizeMosaic();
            end
        end
        
        function reassignMCone(this, ratio, showMosaic)
            if ~exist('showMosaic', 'var')
                showMosaic = true;
            end
            
            this.reassignCone(ratio, this.M_Cone_Idx, this.L_Cone_Idx, showMosaic);
        end
        
        function reassignLCone(this, ratio, showMosaic)
            if ~exist('showMosaic', 'var')
                showMosaic = true;
            end
            
            this.reassignCone(ratio, this.L_Cone_Idx, this.M_Cone_Idx, showMosaic);
        end
        
        function renderMtx = hyperRender(this, imageSize, wave, scale, waitBar)
            if ~exist('waitBar', 'var')
                waitBar = true;
            end
            
            testInput = rand(imageSize);
            testCone = this.computeHyperspectral(wave, scale, testInput);
            
            renderMtx = zeros(length(testCone), length(testInput(:)), 'single');
            updateWaitbar = [];
            if waitBar
                updateWaitbar = waitbarParfor(length(testInput(:)), "Calculation in progress...");
            end
            
            parfor idx = 1:length(testInput(:))
                input = zeros(size(testInput));
                input(idx) = 1.0;
                
                coneVec = this.computeHyperspectral(wave, scale, input);
                renderMtx(:, idx) = single(coneVec);
                
                if waitBar
                    updateWaitbar();
                end
            end
        end
        
        function renderMtx = forwardRender(this, imageSize, validation, optics, waitBar)
            if ~exist('validation', 'var')
                validation = true;
            end
            
            if ~exist('optics', 'var')
                optics = true;
            end
            
            if ~exist('waitBar', 'var')
                waitBar = true;
            end
            
            testInput = rand(imageSize);
            if optics
                [~, ~, testLinear, testCone] = this.compute(testInput);
            else
                [~, ~, testLinear, testCone] = this.computeDiffLmt(testInput);
            end
            
            renderMtx = zeros(length(testCone), length(testLinear(:)), 'single');
            
            updateWaitbar = [];
            if waitBar
                updateWaitbar = waitbarParfor(length(testLinear(:)), "Calculation in progress...");
            end
            parfor idx = 1:length(testLinear(:))
                input = zeros(size(testLinear));
                input(idx) = 1.0;
                
                if optics
                    [~, ~, ~, coneVec] = this.compute(input);
                else
                    [~, ~, ~, coneVec] = this.computeDiffLmt(input);
                end
                
                renderMtx(:, idx) = single(coneVec);
                
                if waitBar
                    updateWaitbar();
                end
            end
            
            if validation
                figure(); hold on;
                scatter(testCone, renderMtx * testLinear(:));
                plot(xlim, ylim, '--k', 'LineWidth', 2);
                axis square;
            end
        end
        
        function reconValidation(this, input, recon, imageSize, coneVec, estimator)
            
            figure();
            plotAxis = tight_subplot(2, 3, [.05 .05], [.05 .05], [.05 .05]);
            
            axes(plotAxis(1));
            imshow(input);
            
            [~, ~, linear, ~] = this.compute(input);
            nlogll = estimator.prior(linear);
            title(sprintf('Original Prior: %.2f \n Loss: %.2f', nlogll, estimator.reconObjective(coneVec, linear(:))));
            axes(plotAxis(2));
            
            oiImage = this.rgbOpticalImage();
            coor = (size(oiImage, 1) - imageSize)/2;
            imshow(imcrop(oiImage, [coor + 1, coor + 1, imageSize, imageSize]));
            title('Optical Image');
            
            axes(plotAxis(3));
            this.visualizeExcitation(true);
            title(sprintf('Cone Excitation \n Likelihood %.2f', estimator.likelihood(coneVec, linear(:))), 'FontSize', 11);
            
            axes(plotAxis(4));
            rgbRecon = invGammaCorrection(recon, this.Display);
            imshow(rgbRecon);
            title('Reconstruction');
            
            [~, ~, linear, ~] = this.compute(rgbRecon);
            nlogll = estimator.prior(linear);
            title(sprintf('Reconstruction: %.2f \n Loss: %.2f', nlogll, estimator.reconObjective(coneVec, linear(:))));
            
            axes(plotAxis(5));
            
            oiImage = this.rgbOpticalImage();
            coor = (size(oiImage, 1) - imageSize)/2;
            imshow(imcrop(oiImage, [coor + 1, coor + 1, imageSize, imageSize]));
            title('Optical Image');
            
            axes(plotAxis(6));
            this.visualizeExcitation(true);
            title(sprintf('Cone Excitation \n Likelihood %.2f', estimator.likelihood(coneVec, linear(:))), 'FontSize', 11);
            
        end
    end
    
    methods (Access = private)
        function response = getConetypeResponse(this, type)
            sizeExci = size(this.LastResponse);
            excitation = reshape(this.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            response   = excitation(this.Mosaic.pattern == type);
        end
    end
    
end

