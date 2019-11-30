classdef PatchEstimator < handle
    
    properties
        Render;  % Render matrix estimated from isetbio rountine
        Basis;   % Image Basis (i.e., PCA or ICA basis)
        Mu;      % Vector of the mean value of each pixel
        Disp;    % Display option for optimization
        Lambda;  % Regularization
        Size;    % Size of original image
        Patch;   % Size of reconstruction blocks
        Stride;  % Step stride size
    end
    
    methods (Abstract)
        
        % Prior on one small patch of the image
        [nlogll, gradient] = priorPatch(this, patchVec)
        
    end
    
    methods
        
        % Constructor for the estimator
        function this = PatchEstimator(render, basis, mu, lambda, stride, imageSize)
            this.Render = render;
            this.Basis  = basis;
            this.Mu     = mu;
            this.Lambda = lambda;
            this.Stride = stride;
            this.Size   = imageSize;
            this.Patch  = sqrt(size(basis, 1) / 3); % Assume square basis image
            this.Disp   = 'iter';
        end
        
        % Gaussian approximation of the likelihood
        function [nlogll, gradient] = likelihood(this, measure, imageVec)
            pred   = this.Render * imageVec;
            nlogll = sum((pred - measure) .^ 2);
            gradient = (2 * (pred - measure)' * this.Render)';
        end
        
        % Prior loss and gradient
        function [nlogll, gradient] = prior(this, image)
            nlogll   = 0;
            gradient = zeros(length(image(:)), 1);
            
            for x = 1:this.Stride:(this.Size(1) - this.Patch + 1)
                for y = 1:this.Stride:(this.Size(2) - this.Patch + 1)
                    idxX = x:1:(x+this.Patch-1);
                    idxY = y:1:(y+this.Patch-1);
                    imagePatch = image(idxX, idxY, :);
                    [nlogllPatch, gradientPatch] = this.priorPatch(imagePatch(:));
                    
                    nlogll = nlogll + nlogllPatch;
                    
                    gradImage = zeros(this.Size);
                    gradImage(idxX, idxY, :) = reshape(gradientPatch, [this.Patch, this.Patch, 3]);
                    
                    gradient = gradient + gradImage(:);
                end
            end
            
            nlogll = this.Lambda * nlogll;
            gradient = this.Lambda * gradient;
        end
        
        function [loss, gradient] = reconObjective(this, measure, imageVec)
            [nlogllPrior, gradientPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogllLlhd,  gradientLlhd]  = this.likelihood(measure, imageVec);
            
            loss = nlogllPrior + nlogllLlhd;
            gradient = gradientPrior + gradientLlhd;
        end
        
        function reconstruction = estimate(this, measure, maxIter, init, bounded, ub)
            loss = @(x) this.reconObjective(measure, x);
            
            if ~exist('maxIter', 'var')
                maxIter = 1e3;
            end
            
            if ~exist('init', 'var')
                init = rand([prod(this.Size), 1]);
            end
            
            if ~exist('bounded', 'var')
                bounded = false;
            end
            
            if ~exist('ub', 'var')
                ub = 1;
            end
            
            if bounded
                options = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxIterations', maxIter, 'CheckGradients', false, ...
                    'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, 'HessianApproximation', 'lbfgs');
                lb = init * 0;
                ub = ones(size(init)) * ub;
                solution = fmincon(loss, init, [], [], [], [], lb, ub, [], options);
            else
                options  = optimset('GradObj', 'on', 'Display', 'iter', 'MaxIter', maxIter);
                solution = fminlbfgs(loss, init, options);
            end
                        
            reconstruction = reshape(solution, this.Size);
        end
        
        function [nlogPrior, gradPrior, nlogLlhd, gradLlhd] = evalEstimate(this, measure, imageVec)
            [nlogPrior, gradPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogLlhd,  gradLlhd]  = this.likelihood(measure, imageVec);
        end
    end
end

