classdef RegressionEstimator < Estimator
%REGRESSIONESTIMATOR Simple regression method for reconstruction.
%
% Syntax: estimator = RegressionEstimator(input, output, [varargin]);
%
% This class implements the linear regression estimator described in Golden
% et al., 2019, with # of diagonal element in SVD as regularization.
%
% RegressionEstimator Methods:
%  estimate   - Compute the linear estimate (e.g., image) based on input
%               (e.g., cone response).
%  setRegPara - Set the regularization parameter.
%
% Inputs:
%   X         - Input (e.g., cone response vector) of the training set
%   Y         - Output (e.g., RGB image) of the training set
%
% Outputs:
%   RegressionEstimator object
%
% Optional key/value pairs:
%   nDiag     - Int. Number of diagonal element used in estimating the weight
%               matrix from SVD.
%               Default is min(szie(X)).
%   paraList  - List of hyperparameter (i.e., number of diagonal elements)
%               we want to test with cross validation.
%               Default is [].   
    
    properties       
        W; % Weight matrix for regression
        X; % Training input cone excitations
        Y; % Training input ground truth images        

        U; 
        invD;
        V;
        
        zeroThreshold = 1e-6;
    end
    
    methods
       
        function obj = RegressionEstimator(X, Y, varargin)
            obj.X = X;
            obj.Y = Y;
            maxDiag = min(size(X));
            
            p = inputParser;
            p.addParameter('nDiag', maxDiag, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('paraList', [], @(x) isnumeric);
            
            parse(p, varargin{:});
            obj.regParaList = p.Results.paraList;
            
            obj.estimateW(p.Results.nDiag);
        end
        
        function reconImage = estimate(obj, input)
            reconImage = input * obj.W;
        end
        
        function estimateW(obj, regPara)
            fprintf('Compute SVD of the input ... ');
            [obj.U, D, obj.V] = svd(obj.X, 'econ');
            fprintf('Done! \n');
            
            diagVec = diag(D);
            nCount  = sum(diagVec <= obj.zeroThreshold);
            
            msg = sprintf('At least %d out of the %d diagonal components are close to zero, \n', ...
                nCount, length(diagVec));
            msg = strcat(msg, ' consider drop them for regularization purpose.');
            
            if(nCount > 0)
                warning(msg);
            end
            
            obj.invD = inv(D);
            obj.setRegPara(regPara);
        end
        
        function setRegPara(obj, regPara)                        
            diagVec = diag(obj.invD);
            diagVec(regPara + 1 : end) = 0;
            Dk = diag(diagVec);
            
            fprintf('Estimate weight matrix ... ');
            obj.W = obj.V * Dk * obj.U' * obj.Y;
            fprintf('Done! \n');
        end
        
    end
    
end