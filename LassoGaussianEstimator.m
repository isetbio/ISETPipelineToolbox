classdef LassoGaussianEstimator < BayesianEstimator
    %LASSOGAUSSIANESTIMATOR
    methods
        function this = LassoGaussianEstimator(render, basis, mu)
            this@BayesianEstimator(render, basis, mu);            
        end
        
        function reconImage = estimate(this, input, varargin)
            p = inputParser;
            p.addParameter('regularization', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            parse(p, varargin{:});
            
            lambda = p.Results.regularization / length(input);
            target = input' - this.combinedBias;
            
            Mdl = fitrlinear(this.combinedRender, target, 'Lambda', lambda, ...
                'Learner', 'leastsquares', 'Regularization', 'lasso', 'FitBias', false, 'Solver', 'sparsa', 'Verbose', 0, ...
                'BetaTolerance', 1e-8, 'GradientTolerance', 1e-8, 'IterationLimit', 5e3);
            reconImage = (this.Basis(:, 1:this.nDim) * Mdl.Beta + this.Mu)';
        end
    end
end

