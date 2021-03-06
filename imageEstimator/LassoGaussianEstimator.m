classdef LassoGaussianEstimator < BayesianEstimator
    %LASSOGAUSSIANESTIMATOR
    properties
        Verbose;
        Lambda;
        IterationLimit;
        Solver;
        Tolerance;
    end
        
    methods
        function this = LassoGaussianEstimator(render, basis, mu)
            this@BayesianEstimator(render, basis, mu);
            this.Lambda = 1;
            this.Solver = 'sparsa';
            this.Tolerance = 1e-6;
            this.IterationLimit = 5e3;            
        end
        
        function reconImage = estimate(this, input, varargin)
            p = inputParser;
            p.addParameter('regularization', this.Lambda, @(x)(isnumeric(x) && numel(x) == 1));
            parse(p, varargin{:});
            
            lambda = p.Results.regularization / length(input);
            target = input' - this.combinedBias;
            
            if(strcmp(this.Solver, 'sparsa'))
                [Mdl, fitInfo] = fitrlinear(this.combinedRender, target, 'Lambda', lambda, ...
                    'Learner', 'leastsquares', 'Regularization', 'lasso', 'FitBias', false, 'Solver', 'sparsa', 'Verbose', ...
                    this.Verbose, 'BetaTolerance', this.Tolerance, 'GradientTolerance', this.Tolerance, 'IterationLimit', this.IterationLimit);
            elseif(strcmp(this.Solver, 'sgd'))
                [Mdl, fitInfo] = fitrlinear(this.combinedRender, target, 'Lambda', lambda, ...
                    'Learner', 'leastsquares', 'Regularization', 'lasso', 'FitBias', false, 'Solver', 'sgd', 'Verbose', ...
                    this.Verbose, 'IterationLimit', this.IterationLimit * length(input), 'BatchSize', length(input), ...
                    'BetaTolerance', 1e-6, 'PassLimit', this.IterationLimit, 'NumCheckConvergence', 5e2);
            else
                error('Invalid Solver Name.')
            end
            fprintf('Search Terminated: %s \n', fitInfo.TerminationStatus{1});
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
        
        function setIterationLimit(this, limit)
            this.IterationLimit = limit;
        end
        
        function setSolver(this, solver)
            this.Solver = solver;
        end
        
        function setTolerance(this, tol)
            this.Tolerance = tol;
        end
    end
end

