classdef BayesianEstimator < Estimator
    
    properties
        Render;  % Render matrix estimated from isetbio rountine
        Basis;   % Image Basis (i.e., PCA or ICA basis)
        Mu;
        
        nDim;
        combinedRender;
        combinedBias;
    end
    
    methods
        
        function obj = BayesianEstimator(render, basis, mu)
            obj.Render = render;
            obj.Basis  = basis;
            obj.Mu     = mu;
            obj.nDim   = size(basis, 1);            
            
            obj.combinedRender = obj.Render * obj.Basis;
            obj.combinedBias   = obj.Render * obj.Mu;            
        end
        
    end
end

