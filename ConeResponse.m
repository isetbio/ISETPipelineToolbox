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
%   visualizeCone        - Visualize the cone mosaic being simulated    
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
                
        function obj = ConeResponse(varargin)
        % CONERESPONSE  Construct ConeResponse object.
        %   Construct ConeResponse object with optional argument:
        %   FovealDegree, eccBasedConeDensity, eccBasedConeQuantal.
            p = inputParser;
            p.addParameter('fovealDegree', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('eccBasedConeDensity', False, @islogical);
            p.addParameter('eccBasedConeQuantal', False, @islogical);
            
            parse(p, varargin);            
            obj.FovealDegree    = p.Results.fovealDegree;
            eccBasedConeDensity = p.Results.eccBasedConeDensity;
            eccBasedConeQuantal = p.Results.eccBasedConeQuantal;                        
            
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
        
        function [excitation, theOI, allCone, L, M, S] = compute(obj, image)
            % COMPUTE   Compute optical image and cone mosaic excitation.                           
            %
            % Syntax: 
            %   [excitation, OI, allCone, L, M, S] = obj.compute(image)
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
            
            allCone = [L; M; S];
                        
        end
        
        function mosaic = getMosaic(obj)
        % Return the cone mosaic object
            mosaic = obj.Mosaic;
        end
        
        function visualizeExcitation(obj)
        % Visualize cone mosaic excitation
            visualizeConeMosaicResponses(obj.Mosaic, obj.LastResponse, 'R*/cone/tau');
        end
        
        function visualizeOI(obj)
        % Visualize optical image
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

