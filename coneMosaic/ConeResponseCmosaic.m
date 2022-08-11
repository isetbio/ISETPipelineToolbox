classdef ConeResponseCmosaic < ConeResponse
    
    properties (Access = public)
        psfData;
        eccX;
        eccY;
    end
    
    methods (Access = public)
        function this = ConeResponseCmosaic(eccX, eccY, varargin)
            p = inputParser;
            p.addParameter('fovealDegree', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('pupilSize', 3.0, @(x) (isnumeric(x) && numel(x) == 1));
            p.addParameter('subjectID', 6, @(x) (isnumeric(x) && numel(x) == 1));
            p.addParameter('randomMesh', false, @islogical);
            p.addParameter('useRandomSeed', true, @islogical)
            
            this@ConeResponse(varargin{:}, 'override', true);
            
            parse(p, varargin{:});
            [mosaic, psfObj, psfData] = PeripheralModel.eyeModelCmosaic...
                (eccX, eccY, p.Results.fovealDegree, p.Results.pupilSize, p.Results.randomMesh, p.Results.subjectID);
            
            this.eccX = eccX;
            this.eccY = eccY;
            
            this.Mosaic = mosaic;
            this.PupilSize = p.Results.pupilSize;
            this.PSF = psfObj;
            this.psfData = psfData;
        end
        
        % Override the visualization with new method
        function visualizeMosaic(this, figureHandle, axesHandle)
            if ~exist('figureHandle', 'var')
                figureHandle = figure();
            end

            if ~exist('axesHandle', 'var')
                axesHandle = [];
            end

            this.Mosaic.visualize('figureHandle', figureHandle, ...
                                  'axesHandle', axesHandle);
        end
        
        function visualizeExcitation(this, figureHandle, axesHandle)
            if ~exist('figureHandle', 'var')
                figureHandle = figure();
            end

            if ~exist('axesHandle', 'var')
                axesHandle = [];
            end

            activationRange = prctile(this.LastResponse(:), [1 99]);
            this.Mosaic.visualize(...                
                'figureHandle', figureHandle, ...
                'axesHandle', axesHandle, ...
                'activation', this.LastResponse, ...
                'activationRange', activationRange, ...
                'plotTitle',  'Cone Response');
        end
        
        function visualizePSF(this)
            [~, wIdx] = min(abs(this.psfData.supportWavelength-550));
            wavePSF = squeeze(this.psfData.data(:,:,wIdx));
            zLevels = 0.1:0.1:0.9;
            xyRangeArcMin = 8 * [-1, 1];
            
            PolansOptics.renderPSF(gca(), ...
                this.psfData.supportX, this.psfData.supportY, wavePSF/max(wavePSF(:)), ...
                xyRangeArcMin, zLevels,  gray(1024), [0 0 0]);
        end
        
        % Compute cone response to image
        function [allCone, linearImage] = compute(this, image)
            % create a visual scene from a display
            meanLuminanceCdPerM2 = [];
            [realizedStimulusScene, ~, linearImage] = sceneFromFile(image, 'rgb', ...
                meanLuminanceCdPerM2, this.Display);
            
            % set the angular scene width
            realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', this.FovealDegree);
            
            % optics
            theOI = oiCompute(realizedStimulusScene, this.PSF);
            this.LastOI = theOI;
            
            % compute cone response
            this.LastResponse = this.Mosaic.compute(theOI, 'opticalImagePositionDegs', 'mosaic-centered');
            allCone = this.LastResponse(:);
        end
        
        function allCone = computeWithScene(this, inputScene)
            opticalImage = oiCompute(inputScene, this.PSF);
            
            % compute cone response
            this.LastOI = opticalImage;
            this.LastResponse = this.Mosaic.compute(opticalImage, 'opticalImagePositionDegs', 'mosaic-centered');
            allCone = this.LastResponse(:);
        end
        
        % Compute render matrix
        function renderMtx = forwardRender(this, imageSize, validation, waitBar, varargin)
            p = inputParser;
            p.addParameter('useDoublePrecision', false, @islogical);
            parse(p, varargin{:});

            if ~exist('validation', 'var')
                validation = true;
            end
            if ~exist('waitBar', 'var')
                waitBar = true;
            end
            
            testInput = rand(imageSize);
            [testCone, testLinear] = this.compute(testInput);
            
            renderMtx = zeros(length(testCone), length(testLinear(:)), 'single');
            
            updateWaitbar = [];
            if waitBar
                updateWaitbar = waitbarParfor(length(testLinear(:)), "Calculation in progress...");
            end
            
            parfor idx = 1:length(testLinear(:))
                input = zeros(size(testLinear));
                input(idx) = 1.0;
                
                coneVec = this.compute(input);
                if (p.Results.useDoublePrecision)
                    renderMtx(:, idx) = single(coneVec);
                else
                    renderMtx(:, idx) = coneVec;
                end
                
                if waitBar
                    updateWaitbar();
                end
            end
            
            if validation
                figure(); hold on;
                scatter(testCone, renderMtx * testLinear(:));
                plot(xlim, xlim, '--k', 'LineWidth', 2);
                axis square;
            end
        end
        
    end
end

