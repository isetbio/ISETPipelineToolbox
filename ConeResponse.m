classdef ConeResponse < handle
    %CONERESPONSE Compute cone responses from RGB image.
    %   This class computes the mean cone excitations and the optical image
    %   for a given mosaic and RGB image. 
    
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
        
        % fovdeg, density, quantal, 
        function obj = ConeResponse(varargin)
            if nargin ~= 0
                obj.FovealDegree    = varargin{1};
                eccBasedConeDensity = varargin{2};
                eccBasedConeQuantal = varargin{3};
            else
                obj.FovealDegree    = 1;
                eccBasedConeDensity = false;
                eccBasedConeQuantal = false;
            end
            
            % Create cone moasic
            fprintf('Create cone moasic object: \n');
            
            theMosaic = coneMosaicHex(5, ...                               % hex lattice sampling factor
                'fovDegs', obj.FovealDegree, ...                           % match mosaic width to stimulus size
                'eccBasedConeDensity', eccBasedConeDensity, ...            % cone density varies with eccentricity
                'eccBasedConeQuantalEfficiency', eccBasedConeQuantal, ...  % cone quantal efficiency varies with eccentricity
                'integrationTime', 0.1, ...                                % 0.1s integration time
                'maxGridAdjustmentIterations', 50);                        % terminate iterative lattice adjustment after 50 iterations            

            % Poisson noise model, mean response
            theMosaic.noiseFlag = 'none';
            obj.Mosaic = theMosaic;

            % Display & Point spread function of human eye
            obj.Display = displayCreate('LCD-Apple', 'viewing distance', 0.57);                       
            obj.PSF = oiCreate('wvf human');            
        end
        
        function visualizeCone(obj)
            obj.Mosaic.visualizeGrid(...
            'backgroundColor', [1 1 1], ...
            'ticksInVisualDegs', true);
        end
        
        function [excitation, theOI, L, M, S] = compute(obj, image)
            % load image
            meanLuminanceCdPerM2 = 100;           
            realizedStimulusScene = sceneFromFile(image, 'rgb', ...
            meanLuminanceCdPerM2, obj.Display);
            
            % set the angular scene width
            realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', obj.FovealDegree);
            
            % optics
            theOI = oiCompute(obj.PSF, realizedStimulusScene);
            obj.LastOI = theOI;
            
            % cone excitations
            nTrialsNum = 2;
            emPath = zeros(nTrialsNum, 1, 2);

            % compute mosaic excitation responses
            % without the eye movement path
            obj.LastResponse = obj.Mosaic.compute(theOI, 'emPath', emPath);             
            
            sizeExci   = size(obj.LastResponse);
            excitation = reshape(obj.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            
            % LMS cone excitation
            L = obj.getConetypeResponse(obj.L_Cone_Idx);
            M = obj.getConetypeResponse(obj.M_Cone_Idx);
            S = obj.getConetypeResponse(obj.S_Cone_Idx);
                        
        end
        
        function visualizeExcitation(obj)
            visualizeConeMosaicResponses(obj.Mosaic, obj.LastResponse, 'R*/cone/tau');
        end
        
        function visualizeOI(obj)
            visualizeOpticalImage(obj.LastOI, 'displayRetinalContrastProfiles', true);
        end
                         
    end    
    
    methods (Access = private)
        
        function response = getConetypeResponse(obj, type)
            sizeExci = size(obj.LastResponse);
            excitation = reshape(obj.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
            response   = excitation(obj.Mosaic.pattern == type);
        end
    end
    
end

