classdef LassoDenoiseEstimator < BayesianEstimator
    %LASSOESTIMATOR
    properties
        invRender;
    end
    
    methods
        function this = LassoDenoiseEstimator(render, basis, mu)
            this@BayesianEstimator(render, basis, mu);
            this.invRender = inv(render);
        end
        
        function reconImage = estimate(this, input, varargin)
            p = inputParser;
            p.addParameter('regularization', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            parse(p, varargin{:});
                        
            alpha  = 1;
            lambda = p.Results.regularization / (2 * length(input));
            target = input' - this.combinedBias;
            
            [coff, fitInfo] = lasso(this.combinedRender, target, ...
                'Alpha', alpha, 'Lambda', lambda, 'MaxIter', 1e6);
            reconImage = (this.Basis(:, 1:this.nDim) * coff + this.Mu + ...
                this.invRender * fitInfo.Intercept * ones(length(input), 1))';
        end
    end
end

