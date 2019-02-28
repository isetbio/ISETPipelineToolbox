classdef PoissonGaussianEst < Estimator
    %POISSONGAUSSIANEST
    
    properties        
        Render;     % Render matrix estimated from isetbio rountine
        Basis;      % PCA basis
        Mu;         % PCA mu vector
        
        nDim;
        combinedRender;
        combinedBias;
    end
    
    methods
        
        function obj = PoissonGaussianEst(render, basis, mu)
            obj.Render = render;
            obj.Basis  = basis;
            obj.Mu     = mu;
            obj.nDim   = size(basis, 1);
            
            obj.combinedRender = obj.Render * obj.Basis;
            obj.combinedBias   = obj.Render * obj.Mu;
        end
        
        % Poisson log likelihood
        function logll = logll(~, excitation, lambda)                     
            idpdLl = -lambda + excitation' .* log(lambda);            
            logll = sum(idpdLl);
        end
        
        % Gradient of likelihood
        function gradll = gradll(obj, excitation, x)
            lambda = obj.combinedRender * x + obj.combinedBias;
            dldx1  = obj.combinedRender;
            dldx2  = (-excitation') .* obj.combinedRender ./ lambda;
            
            gradll = sum(dldx1 + dldx2, 1);
        end
        
        % Gradient of prior
        function gradprior = gradprior(~, x)            
        end
        
        % Posterior likelihood
        function negll = negll(obj, input, x)
            priorLoss = - sum(log(normpdf(x, 0, 1)));
            logllLoss = - obj.logll(input, obj.combinedRender * x + obj.combinedBias);
            negll = priorLoss + logllLoss;
        end
        
        function reconImage = estimate(obj, input)
            loss = @(x) obj.negll(input, x);
                                                
            % Optimization
            problem = createOptimProblem('fmincon');
            problem.objective = loss;
            problem.x0 = zeros([obj.nDim, 1]);
            problem.A = [obj.Basis(:, 1:obj.nDim), -obj.Basis(:, 1:obj.nDim)];
            problem.b = [1 - obj.Mu; obj.Mu];
            problem.Display = 'iter';            
            
            coff = fmincon(problem);            
            reconImage = obj.combinedRender * coff + obj.combinedBias;
        end
        
        function setRegPara(obj, para)
            obj.nDim = para;
            obj.combinedRender = obj.Render * obj.Basis(:, 1:obj.nDim);
        end
        
    end
    
end

