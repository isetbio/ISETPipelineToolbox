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
            assert(sum(lambda <= 0) == 0, 'lambda <= 0');
            idpdLl = -lambda + excitation' .* log(lambda);
            logll = sum(idpdLl);
        end
        
        % Gradient of likelihood
        function gradll = gradll(obj, excitation, lambda)                                    
            dldx1  = obj.combinedRender;
            dldx2  = (-excitation') .* obj.combinedRender ./ lambda;
            
            gradll = sum(dldx1 + dldx2, 1);
        end                
        
        % Posterior likelihood and gradient of posterir likelihood
        function [negll, grad] = negll(obj, input, x)
            priorLoss = - sum(log(normpdf(x, 0, 1)));
            
            lambda = obj.combinedRender * x + obj.combinedBias;
            logllLoss = - obj.logll(input, lambda);
            negll = priorLoss + logllLoss;
            
            grad = obj.gradll(input, lambda) + x;
        end
        
        function reconImage = estimate(obj, input)
            loss = @(x) obj.negll(input, x);
                                                
            % Optimization
            problem = createOptimProblem('fmincon');
            problem.objective = loss;
            problem.x0 = zeros([obj.nDim, 1]);
            problem.Aineq = [obj.Basis(:, 1:obj.nDim), -obj.Basis(:, 1:obj.nDim)];
            problem.bineq = [1 - obj.Mu; obj.Mu];                         
            problem.options = ...
                optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1e6, ...
                'SpecifyObjectiveGradient', true, 'CheckGradients',true);            
            
            coff = fmincon(problem);            
            reconImage = obj.Basis(:, 1:obj.nDim) * coff + obj.Mu;
        end
        
        function setRegPara(obj, para)
            obj.nDim = para;
            obj.combinedRender = obj.Render * obj.Basis(:, 1:obj.nDim);
        end
        
    end
    
end

