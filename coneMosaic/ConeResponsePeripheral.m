classdef ConeResponsePeripheral < ConeResponse
    
    methods (Access = public)
        function this = ConeResponsePeripheral(eccX, eccY, varargin)
            display = displayCreate('LCD-Apple');
            
            p = inputParser;
            p.addParameter('fovealDegree', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('pupilSize', 3.0, @(x) (isnumeric(x) && numel(x) == 1));
            p.addParameter('display', display);
            
            this@ConeResponse(varargin{:}, 'override', true);
            
            parse(p, varargin{:});
            [mosaic, psf] = PeripheralModel.eyeModelEcc(eccX, eccY, p.Results.fovealDegree, p.Results.pupilSize);
            
            mosaic.noiseFlag = 'none';
            mosaic.integrationTime = 0.2;
            this.Mosaic = mosaic;
            this.DefaultMosaic = this.Mosaic.pattern;
            
            this.PupilSize = p.Results.pupilSize;
            this.PSF = psf;
        end
    end
end

