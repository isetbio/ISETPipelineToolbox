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
            p.addParameter('useRandomSeed', true, @islogical);
            p.addParameter('defocusDiopters', 0, @isnumeric);
            p.addParameter('wave',400:10:700, @isnumeric);
            p.addParameter('pigment', cPhotoPigment(), @(x) isa(x, 'cPhotoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('tritanopicRadiusDegs', 0.15, @isscalar);
            p.addParameter('rodIntrusionAdjustedConeAperture', false, @(x) ((islogical(x))||((isscalar(x))&&((x>0)&&(x<=1)))));
            p.addParameter('eccVaryingConeAperture', true, @islogical);
            p.addParameter('eccVaryingConeBlur', false, @islogical);
            p.addParameter('eccVaryingOuterSegmentLength', true, @islogical);
            p.addParameter('eccVaryingMacularPigmentDensity', true, @islogical);
            p.addParameter('eccVaryingMacularPigmentDensityDynamic', false, @islogical);
            p.addParameter('anchorAllEccVaryingParamsToTheirFovealValues', false, @islogical);
            p.addParameter('zernikeDataBase', 'Polans2015', @ischar);
            p.addParameter('noLCA',false,@islogical);
            parse(p, varargin{:});

            % Setting override to true here avoids older coneMosaicHex
            % code, and allows us to set the cMosaic below.
            this@ConeResponse(varargin{:}, 'override', true);
            
            [mosaic, psfObj, psfData] = PeripheralModel.eyeModelCmosaic ...
                (eccX, eccY, p.Results.fovealDegree, p.Results.pupilSize, p.Results.randomMesh, p.Results.subjectID, ...
                'useRandomSeed',p.Results.useRandomSeed,'defocusDiopters',p.Results.defocusDiopters,'wave',p.Results.wave, ...
                'tritanopicRadiusDegs', p.Results.tritanopicRadiusDegs, ...
                'macular', p.Results.macular, ...      
                'pigment', p.Results.pigment, ...
                'rodIntrusionAdjustedConeAperture', p.Results.rodIntrusionAdjustedConeAperture, ...
                'eccVaryingConeAperture', p.Results.eccVaryingConeAperture, ...
                'eccVaryingConeBlur', p.Results.eccVaryingConeBlur, ...
                'eccVaryingOuterSegmentLength', p.Results.eccVaryingOuterSegmentLength, ...
                'eccVaryingMacularPigmentDensity', p.Results.eccVaryingMacularPigmentDensity, ...
                'eccVaryingMacularPigmentDensityDynamic', p.Results.eccVaryingMacularPigmentDensity, ...
                'anchorAllEccVaryingParamsToTheirFovealValues', p.Results.anchorAllEccVaryingParamsToTheirFovealValues, ...
                'zernikeDataBase',p.Results.zernikeDataBase, ...
                'noLCA',p.Results.noLCA);
            
            this.eccX = eccX;
            this.eccY = eccY;
            
            this.Mosaic = mosaic;
            this.PupilSize = p.Results.pupilSize;
            this.PSF = psfObj;
            this.psfData = psfData;
        end
        
        % Override the visualization with new method
        function visualizeMosaic(this, figureHandle, axesHandle, domainVisualizationLimits)
            if ~exist('figureHandle', 'var')
                figureHandle = figure();
            end

            if ~exist('axesHandle', 'var')
                axesHandle = [];
            end

            if ~exist('domainVisualizationLimits', 'var')
                domainVisualizationLimits = [];
            end

            this.Mosaic.visualize('figureHandle', figureHandle, ...
                                  'axesHandle', axesHandle, ...
                                  'domainVisualizationLimits', domainVisualizationLimits);
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
        function [allCone, linearStimulusImage] = compute(this, stimulusImageRGB, varargin)

            p = inputParser;
            p.addParameter('lowOpticalImageResolutionWarning', true, @islogical);
            p.parse(varargin{:});

            % Create a visual scene from a display
            meanLuminanceCdPerM2 = [];
            [stimulusScene, ~, linearStimulusImage] = sceneFromFile(stimulusImageRGB, 'rgb', ...
                meanLuminanceCdPerM2, this.Display);
            
            % Set the angular scene width
            stimulusScene = sceneSet(stimulusScene, 'fov', this.FovealDegree);
            
            % optics
            theOI = oiCompute(stimulusScene, this.PSF);
            this.LastOI = theOI;
            
            % compute cone response
            this.LastResponse = this.Mosaic.compute(theOI, 'opticalImagePositionDegs', 'mosaic-centered', ...
                'lowOpticalImageResolutionWarning',p.Results.lowOpticalImageResolutionWarning);
            allCone = this.LastResponse(:);
        end
        
        function allCone = computeWithScene(this, inputScene)
            opticalImage = oiCompute(inputScene, this.PSF);
            
            % compute cone response
            this.LastOI = opticalImage;
            this.LastResponse = this.Mosaic.compute(opticalImage, 'opticalImagePositionDegs', 'mosaic-centered');
            allCone = this.LastResponse(:);
        end
        
        % Compute render matrix.  Each monitor channel is sequential in the
        % column vecttor the represents an image.
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
            
            if (p.Results.useDoublePrecision)
                renderMtx = zeros(length(testCone), length(testLinear(:)));
            else
                renderMtx = zeros(length(testCone), length(testLinear(:)), 'single');
            end
            
            updateWaitbar = [];
            if waitBar
                updateWaitbar = waitbarParfor(length(testLinear(:)), "Calculation in progress...");
            end
            
            input = zeros(size(testLinear));
                input(1) = 1.0;

            parfor idx = 1:length(testLinear(:))
                input = zeros(size(testLinear));
                input(idx) = 1.0;
                
                % Only throw the optical image resolution warning once if
                % it is thrown for the current configuration.
                [coneVec, chkLinear] = this.compute(input,'lowOpticalImageResolutionWarning',false);
                if (~p.Results.useDoublePrecision)
                    renderMtx(:, idx) = single(coneVec);
                else
                    renderMtx(:, idx) = coneVec;
                end
                if (any(input(:) ~= chkLinear(:)))
                    error('Conversion to linear not as expected for RGB input of 1\n');
                end

                % You can uncomment this line in a desparate attempt to figure
                % out why render matrix doesn't reproduce test responses.
                %
                % if (max(coneVec(:)) < 10)
                %     fprintf('Warning: Small maximum cone excitation %g\n',max(coneVec(:)));
                % end
                
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
        
        % Render matrix version where each rgb channel of monitor
        % representation has been scaled by the corresponding entry
        % of scaleVector. Used to tweak monitor scaling without needing
        % to recompute all of the nice render matrices we may already have.
        function renderMtx = scaleRenderMtx(this, renderMtx, scaleVector, varargin)

            [m,n] = size(renderMtx);
            nPrimaries = length(scaleVector);
            nPixels = n/nPrimaries;
            if (nPixels ~= floor(nPixels))
                error('Logic error recreating number of pixels from render matrix');
            end

            for jj = 1:nPixels
                for pp = 1:nPrimaries
                   renderMtx(:,jj+(pp-1)*nPixels) = renderMtx(:,jj+(pp-1)*nPixels)*scaleFactor(pp)
                end
            end

        end
    end
end

