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
            zeroBnd = 1e-5;
            lambda(lambda < zeroBnd) = zeroBnd;
            
            idpdLl = -lambda + excitation' .* log(lambda);            
            logll = sum(idpdLl);
        end
        
        % Posterior likelihood
        function negll = negll(obj, input, x)
            priorLoss = - sum(log(normpdf(x, 0, 1)));
            logllLoss = - obj.logll(input, obj.combinedRender * x + obj.combinedBias);
            negll = priorLoss + logllLoss;
        end
        
        function reconImage = estimate(obj, input)
            loss = @(x) obj.negll(input, x);
            
            init = normrnd(0, 1, [obj.nDim, 1]);
            % Optimization
            options = optimoptions('fminunc');
            options.MaxFunctionEvaluations = 1e6;
            options.Display = 'iter';            
            coff = fminunc(loss, init, options);
            
            reconImage = obj.combinedRender * coff + obj.combinedBias;
        end
        
        function setRegPara(obj, para)
            obj.nDim = para;
            obj.combinedRender = obj.Render * obj.Basis(:, 1:obj.nDim);
        end
        
    end
    
end

