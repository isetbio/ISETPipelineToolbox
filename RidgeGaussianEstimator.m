classdef RidgeGaussianEstimator < BayesianEstimator
    %RIDGEGAUSSIANESTIMATOR 
       
    methods
        
        function obj = RidgeGaussianEstimator(render, basis, mu)
            obj@BayesianEstimator(render, basis, mu);
        end
        
        function reconImage = estimate(obj, input)
            coff = (obj.combinedRender' * obj.combinedRender + ...
                diag(ones(1, obj.nDim))) \ obj.combinedRender' * (input' - obj.combinedBias);
            reconImage = (obj.Basis(:, 1:obj.nDim) * coff + obj.Mu)';            
        end
        
    end
end

