classdef GaussianPatchEstimator < PatchEstimator
    
    methods
        
        % Constructor for the estimator
        function this = GaussianPatchEstimator(render, basis, mu, lambda, stride, imageSize)
            this@PatchEstimator(render, basis, mu, lambda, stride, imageSize);
        end
        
        % Patchwise Gaussian prior and gradient
        function [nlogll, gradient] = priorPatch(this, patchVec)
            projection = this.Basis * (patchVec - this.Mu);
            nlogll   = 0.5 * sum(projection .^ 2);
            gradient = projection' * this.Basis;
        end
               
    end
    
end

