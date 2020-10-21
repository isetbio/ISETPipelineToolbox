classdef RidgeGaussianEstimator < BayesianEstimator
    %RIDGEGAUSSIANESTIMATOR
    properties
        Lambda;
    end
    
    methods
        
        function this = RidgeGaussianEstimator(render, basis, mu)
            this@BayesianEstimator(render, basis, mu);
            this.Lambda = 1;
        end
        
        function reconImage = estimate(this, input, varargin)
            p = inputParser;
            p.addParameter('regularization', this.Lambda, @(x)(isnumeric(x) && numel(x) == 1));
            parse(p, varargin{:});
            
            coff = (this.combinedRender' * this.combinedRender + ...
                p.Results.regularization * diag(ones(1, this.nDim))) \ this.combinedRender' * (input' - this.combinedBias);
            reconImage = (this.Basis(:, 1:this.nDim) * coff + this.Mu)';
        end
        
        function setLambda(this, lambda)
            this.Lambda = lambda;
        end
        
    end
end

