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
        Mosaic;
        LastResponse;
        LastOI;
    end
    
    properties (Access = private)
        L_Cone_Idx = 2;
        M_Cone_Idx = 3;
        S_Cone_Idx = 4;
    end
    
    methods (Access = public)
        
        function this = ConeResponse(varargin)
            % CONERESPONSE  Construct ConeResponse object.
            %   Construct ConeResponse object with optional argument:
            %   FovealDegree, eccBasedConeDensity, eccBasedConeQuantal.
            display = displayCreate('LCD-Apple');
            
            p = inputParser;
            p.addParameter('fovealDegree', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('eccBasedConeDensity', false, @islogical);
            p.addParameter('eccBasedConeQuantal', false, @islogical);
            p.addParameter('display', display);
            p.addParameter('viewDistance', 0.57, @(x) (isnumeric(x) && numel(x) == 1));
            p.addParameter('spatialDensity', []);
            
            parse(p, varargin{:});
            this.FovealDegree    = p.Results.fovealDegree;
            spatialDensity      = p.Results.spatialDensity;
            eccBasedConeDensity = p.Results.eccBasedConeDensity;
            eccBasedConeQuantal = p.Results.eccBasedConeQuantal;
            
            % Create cone moasic
            fprintf('Create cone moasic object: \n');
            
            if(isempty(spatialDensity))
                theMosaic = coneMosaicHex(5, ...                               % hex lattice sampling factor
                    'fovDegs', this.FovealDegree, ...                           % match mosaic width to stimulus size
                    'eccBasedConeDensity', eccBasedConeDensity, ...            % cone density varies with eccentricity
                    'eccBasedConeQuantalEfficiency', eccBasedConeQuantal, ...  % cone quantal efficiency varies with eccentricity
                    'integrationTime', 0.2, ...                                % 0.1s integration time
                    'maxGridAdjustmentIterations', 50);                        % terminate iterative lattice adjustment after 50 iterations
            else
                theMosaic = coneMosaicHex(5, ...
                    'fovDegs', this.FovealDegree, ...
                    'spatialDensity', spatialDensity, ...
                    'sConeMinDistanceFactor', 0, ...
                    'sConeFreeRadiusMicrons', 0, ...
                    'integrationTime', 0.2, ...
                    'maxGridAdjustmentIterations', 50);
            end
            
            % Poisson noise model, mean response
            theMosaic.noiseFlag = 'none';
            this.Mosaic = theMosaic;
            
            % Display
            this.Display = p.Results.display;
            this.Display.dist = p.Results.viewDistance;
            
            % Point spread function of human eye
            this.PSF = oiCreate('wvf human');
        end
        
        function visualizeMosaic(this)
            this.Mosaic.visualizeGrid(...
                'backgroundColor', [1 1 1], ...
                'ticksInVisualDegs', true);
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
        
        function [allCone] = coneExcitationRnd(this, factor, type)
            sizeExci   = size(this.LastResponse);
            excitation = reshape(this.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            
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
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', 'R*/cone/tau');
            set(gca, 'XTick', [], 'YTick', []);
            set(gca,'YDir','reverse');
            title('Cone Excitation Pattern');
        end
        
        function visualizeOI(this)
            % Visualize optical image
            visualizeOpticalImage(this.LastOI, 'displayRadianceMaps', false, ...
                'displayRetinalContrastProfiles', false);
            set(gca,'YDir','reverse');
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

