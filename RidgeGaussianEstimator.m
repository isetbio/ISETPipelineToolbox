classdef RidgeGaussianEstimator < BayesianEstimator
    %RIDGEGAUSSIANESTIMATOR 
       
    methods
        
        function obj = RidgeGaussianEstimator(render, basis, mu)
            obj@BayesianEstimator(render, basis, mu);
        end
        
        function reconImage = estimate(obj, input, varargin)
            p = inputParser;
            p.addParameter('regularization', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            parse(p, varargin{:});
            
            coff = (obj.combinedRender' * obj.combinedRender + ...
                p.Results.regularization * diag(ones(1, obj.nDim))) \ obj.combinedRender' * (input' - obj.combinedBias);
            reconImage = (obj.Basis(:, 1:obj.nDim) * coff + obj.Mu)';            
        end
        
    end
end

