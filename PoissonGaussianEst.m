classdef PoissonGaussianEst < Estimator
    %POISSONGAUSSIANEST
    
    properties        
        Render;     % Render matrix estimated from isetbio rountine
        Basis;      % PCA basis
        Mu;         % PCA mu vector
        
        nDim;
    end
    
    methods
        
        function obj = PoissonGaussianEst(render, basis, mu)
            obj.Render = render;
            obj.Basis  = basis;
            obj.Mu     = mu;
            obj.nDim   = size(basis, 1);
        end
        
        % Poisson log likelihood
        function logll = logll(obj, excitation, image)
            lambda = obj.Render * image;
            idpdLl = -lambda + excitation .* log(lambda);
            
            logll = sum(idpdLl);
        end
        
        % Posterior likelihood
        function negll = negll(obj, input, x)
            priorLoss = - sum(log(normpdf(x, 0, 1)));
            logllLoss = - obj.logll(input, obj.Basis * x + obj.Mu);
            negll = priorLoss + logllLoss;
        end
        
        function reconImage = estimate(obj, input)
            loss = @(x) obj.negll(input, x);
            
            [~, dim] = size(pcaBasis);
            init = normrnd(0, 1, [dim, 1]);

            % Optimization
            options = optimoptions('fminunc');
            options.MaxFunctionEvaluations = 1e6;
            options.Display = 'iter';
            options.UseParallel = true;
            coff = fminunc(loss, init, options);
            
            reconImage = obj.Basis * coff + obj.Mu;
        end
        
        function setRegPara(obj, para)
            obj.nDim = para;
        end
        
    end
    
end

