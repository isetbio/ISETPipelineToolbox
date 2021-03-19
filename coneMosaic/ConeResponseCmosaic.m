classdef ConeResponseCmosaic < ConeResponse
    
    properties (Access = public)
        psfData;
    end
    
    methods (Access = public)
        function this = ConeResponseCmosaic(eccX, eccY, varargin)
            p = inputParser;
            p.addParameter('fovealDegree', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('pupilSize', 3.0, @(x) (isnumeric(x) && numel(x) == 1));
            
            this@ConeResponse(varargin{:}, 'override', true);
            
            parse(p, varargin{:});
            [mosaic, psfObj, psfData] = PeripheralModel.eyeModelCmosaic(eccX, eccY, p.Results.fovealDegree, p.Results.pupilSize);
            
            this.Mosaic = mosaic;
            this.PupilSize = p.Results.pupilSize;
            this.PSF = psfObj;
            this.psfData = psfData;
        end
        
        % Override the visualization with new method
        function visualizeMosaic(this)
            this.Mosaic.visualize('figureHandle', figure());
        end
        
        function visualizeExcitation(this)
            activationRange = prctile(this.LastResponse(:), [1 99]);
            this.Mosaic.visualize('figureHandle', figure(), ...
                'activation', this.LastResponse, ...
                'activationRange', activationRange, ...                
                'plotTitle',  'Cone Response');
        end
        
        function visualizePSF(this)            
            [~, wIdx] = min(abs(this.psfData.supportWavelength-550));
            wavePSF = squeeze(this.psfData.data(:,:,wIdx));
            zLevels = 0.1:0.1:0.9;
            xyRangeArcMin = 3;
            
            figure();
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
            theOI = oiCompute(this.PSF, realizedStimulusScene);
            this.LastOI = theOI;
            
            % compute cone response
            this.LastResponse = this.Mosaic.compute(theOI, 'nTrials', 1);
            allCone = this.LastResponse(:);
        end
    end
end

