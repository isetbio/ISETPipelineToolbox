classdef RidgeGaussianEstimator < BayesianEstimator
    %RIDGEGAUSSIANESTIMATOR 
       
    methods
        
        function obj = RidgeGaussianEstimator(render, basis, mu)
            obj@BayesianEstimator(render, basis, mu);
        end
        
        function reconImage = estimate(obj, input, reg)
            if ~exist('reg', 'var')
                reg = 1;
            end            
            coff = (obj.combinedRender' * obj.combinedRender + ...
                reg * diag(ones(1, obj.nDim))) \ obj.combinedRender' * (input' - obj.combinedBias);
            reconImage = (obj.Basis(:, 1:obj.nDim) * coff + obj.Mu)';            
        end
        
    end
end

