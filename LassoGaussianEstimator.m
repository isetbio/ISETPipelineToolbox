classdef LassoGaussianEstimator < BayesianEstimator
    %LASSOGAUSSIANESTIMATOR
    properties
        Verbose;
        Lambda;
    end
        
    methods
        function this = LassoGaussianEstimator(render, basis, mu)
            this@BayesianEstimator(render, basis, mu);
            this.Lambda = 1;
        end
        
        function reconImage = estimate(this, input, varargin)
            p = inputParser;
            p.addParameter('regularization', this.Lambda, @(x)(isnumeric(x) && numel(x) == 1));
            parse(p, varargin{:});
            
            lambda = p.Results.regularization / length(input);
            target = input' - this.combinedBias;
            
            Mdl = fitrlinear(this.combinedRender, target, 'Lambda', lambda, ...
                'Learner', 'leastsquares', 'Regularization', 'lasso', 'FitBias', false, 'Solver', 'sparsa', 'Verbose', ...
                this.Verbose, 'BetaTolerance', 1e-8, 'GradientTolerance', 1e-8, 'IterationLimit', 5e3);
            reconImage = (this.Basis(:, 1:this.nDim) * Mdl.Beta + this.Mu)';
        end
        
        function dispOn(this)
            this.Verbose = 1;
        end
        
        function dispOff(this)
            this.Verbose = 0;
        end
        
        function setLambda(this, lambda)
            this.Lambda = lambda;
        end
    end
end

