classdef SparsePatchEstimator < PatchEstimator

    methods
            
        % Constructor for the estimator
        function this = SparsePatchEstimator(render, basis, mu, lambda, stride, imageSize)
            this@PatchEstimator(render, basis, mu, lambda, stride, imageSize);
        end
        
        % Patchwise Gaussian prior and gradient
        function [nlogll, gradient] = priorPatch(this, patchVec)
            projection = this.Basis * (patchVec - this.Mu);
            nlogll   = sum(abs(projection));

            zero = 1e-2;
            gradient = ((projection > zero) - (projection < -zero))' * this.Basis;             
        end
            
    end

end
