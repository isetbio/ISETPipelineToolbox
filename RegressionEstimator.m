classdef RegressionEstimator < Estimator
    %REGRESSIONESTIMATOR Simple regression method for reconstruction    
    
    properties       
        W; % Weight matrix for regression
        X; % Training input cone excitations
        Y; % Training input ground truth images        

        U;
        D;
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
            [obj.U, obj.D, obj.V] = svd(obj.X, 'econ');
            fprintf('Done! \n');
            
            diagVec = diag(obj.D);
            nCount  = sum(diagVec <= obj.zeroThreshold);
            
            msg = sprintf('At least %d out of the %d diagonal component is too small, \n', nCount, length(diagVec));
            msg = strcat(msg, ' consider drop them for regularization purpose.');
            
            if(nCount > 0)
                warning(msg);
            end
            
            obj.setRegPara(regPara);
        end
        
        function setRegPara(obj, regPara)                        
            diagVec = diag(inv(obj.D));
            diagVec(regPara + 1 : end) = 0;
            Dk = diag(diagVec);
            
            fprintf('Estimate weight matrix ... ');
            obj.W = obj.V * Dk * obj.U' * obj.Y;
            fprintf('Done! \n');
        end
        
    end
    
end