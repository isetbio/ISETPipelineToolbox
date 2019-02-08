classdef RegressionEstimator < Estimator
    %REGRESSIONESTIMATOR Simple regression method for reconstruction    
    
    properties       
        W; % Weight matrix for regression
        X; % Training input cone excitations
        Y; % Training input ground truth images
                
        U;
        D;
        V;
    end
    
    methods
       
        function obj = RegressionEstimator(X, Y, varargin)
            obj.X = X;
            obj.Y = Y;
            
        end
        
        function estimateW(obj, regPara)
            [obj.U, obj.D, obj.V] = svd(obj.X, 'econ');
            obj.setRegPara(regPara);                        
        end
        
        function setRegPara(obj, regPara)            
            diagVec = diag(inv(obj.D));
            diagVec(regPara + 1 : end) = 0;
            Dk = diag(diagVec);
            
            obj.W = obj.V * Dk * obj.U' * obj.Y;
        end
        
    end
    
end