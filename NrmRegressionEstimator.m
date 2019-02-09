classdef NrmRegressionEstimator < RegressionEstimator
    %NRMREGRESSIONESTIMATOR Regression method for reconstruction with
    %normalization as preprocessing step.    
    
    properties
        X_mu;
        X_sigma;
        
        Y_mu;
        Y_sigma;
    end
    
    methods
        
        function obj = NrmRegressionEstimator(X, Y, varargin)
            [X_zscored, X_mu, X_sigma] = zscore(X);
            [Y_zscored, Y_mu, Y_sigma] = zscore(Y);
            
            obj@RegressionEstimator(X_zscored, Y_zscored, varargin{:});
            
            obj.X_mu = X_mu; obj.X_sigma = X_sigma;
            obj.Y_mu = Y_mu; obj.Y_sigma = Y_sigma;
        end
        
        function reconImage = estimate(obj, input)
            input = (input - obj.X_mu) ./ obj.X_sigma;
            reconImage = (input * obj.W) .* obj.Y_sigma + obj.Y_mu;            
        end
        
    end
    
end

