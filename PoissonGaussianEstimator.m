classdef PoissonGaussianEstimator < BayesianEstimator
    %POISSONGAUSSIANEST
    
    properties
        Disp;       % Display option for optimization                
    end
    
    methods
        
        function obj = PoissonGaussianEstimator(render, basis, mu)
            obj@BayesianEstimator(render, basis, mu);
            obj.Disp   = 'iter';
        end
        
        % Poisson log likelihood
        function logll = logll(~, excitation, lambda)            
            idpdLl = -lambda + excitation' .* log(lambda);
            logll = sum(idpdLl);
        end
        
        % Gradient of likelihood
        function gradll = gradll(obj, excitation, lambda)                                    
            dldx1  = obj.combinedRender;
            dldx2  = (-excitation') .* obj.combinedRender ./ lambda;
            
            gradll = sum(dldx1 + dldx2, 1);
        end                
        
        % Posterior likelihood and gradient of posterior likelihood
        function [negll, grad] = negll(obj, input, x)
            priorLoss = - sum(log(normpdf(x, 0, 1)));
            
            lambda = obj.combinedRender * x + obj.combinedBias;
            logllLoss = - obj.logll(input, lambda);
            
            negll = priorLoss + logllLoss;            
            grad = obj.gradll(input, lambda)' + x;            
        end
        
        function reconImage = estimate(obj, input)
            loss = @(x) obj.negll(input, x);
                                                
            % Optimization
            problem = createOptimProblem('fminunc');
            problem.objective = loss;
            problem.x0 = zeros([obj.nDim, 1]);                        
            problem.options = ...
                optimoptions('fminunc', 'Display', obj.Disp, 'SpecifyObjectiveGradient', true, ...
                 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 5e2);
            
            coff = fminunc(problem);
            reconImage = (obj.Basis(:, 1:obj.nDim) * coff + obj.Mu)';
        end                
        
        function dispOn(obj)
            obj.Disp = 'iter';
        end
        
        function dispFinal(obj)
            obj.Disp = 'final';
        end
        
        function dispOff(obj)
            obj.Disp = 'off';
        end
        
    end
    
end

