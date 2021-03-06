classdef PoissonSparseEstimator < SparsePatchEstimator
    
    methods
        % Constructor for the estimator
        function this = PoissonSparseEstimator(render, basis, mu, lambda, stride, imageSize)
            this@SparsePatchEstimator(render, basis, mu, lambda, stride, imageSize);
        end
        
        % Override the likelihood function
        function [nlogll, gradient] = likelihood(this, measure, imageVec)
            lambda = this.Render * imageVec;
            idpdLl = -lambda + measure .* log(lambda);
            nlogll = -sum(idpdLl);
            gradient = (sum(this.Render + (-measure) .* this.Render ./ lambda, 1))';
        end
    end
end

