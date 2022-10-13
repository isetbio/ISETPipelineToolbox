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
        LossFactor % Normalizing factor for loss
        LossMultiplier % Value of loss function at initizialization, after normalizing
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
            this.LossFactor = 1;
            this.LossMultiplier = 10;
        end

        % Gaussian approximation of the likelihood
        function [nlogll, gradient] = likelihood(this, measure, imageVec)
            pred   = this.Render * imageVec;
            nlogll = sum((pred - measure) .^ 2);
            gradient = (2 * (pred - measure)' * this.Render)';
        end

        % Prior loss and gradient
        function [nlogprior, gradient] = prior(this, image)
            nlogprior   = 0;
            gradient = zeros(length(image(:)), 1);

            for x = 1:this.Stride:(this.Size(1) - this.Patch + 1)
                for y = 1:this.Stride:(this.Size(2) - this.Patch + 1)
                    idxX = x:1:(x+this.Patch-1);
                    idxY = y:1:(y+this.Patch-1);
                    imagePatch = image(idxX, idxY, :);
                    [nlogllPatch, gradientPatch] = this.priorPatch(imagePatch(:));

                    nlogprior = nlogprior + nlogllPatch;

                    gradImage = zeros(this.Size);
                    gradImage(idxX, idxY, :) = reshape(gradientPatch, [this.Patch, this.Patch, 3]);

                    gradient = gradient + gradImage(:);
                end
            end

            nlogprior = this.Lambda * nlogprior;
            gradient = this.Lambda * gradient;
        end

        function [loss, gradient] = reconObjective(this, measure, imageVec)
            [nlogllPrior, gradientPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogllLlhd,  gradientLlhd]  = this.likelihood(measure, imageVec);

            loss = nlogllPrior + nlogllLlhd;
            gradient = gradientPrior + gradientLlhd;

            % Scale loss
            loss = this.LossFactor*loss;
            gradient = this.LossFactor*gradient;
        end

        % use GPU (gpuArray) to compute the large matrix product in the
        % likelihood function calculation
        function [loss, gradient] = reconObjectiveGPU(this, measure, imageVec)
            [nlogllPrior, gradientPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogllLlhd,  gradientLlhd]  = this.likelihood(measure, gpuArray(single(imageVec)));

            loss = nlogllPrior + double(gather(nlogllLlhd));
            gradient = gradientPrior + double(gather(gradientLlhd));

            % Scale loss
            loss = this.LossFactor*loss;
            gradient = this.LossFactor*gradient;
        end

        % Run estimate from multiple starting points.
        % 
        % 
        function [multistartStruct,reconImageLinear,reconIndex] = runMultistartEstimate(this, coneVec, varargin)
            p = inputParser;
            p.KeepUnmatched = false;
            p.addParameter('nWhiteStart', 0, @isnumeric);
            p.addParameter('nPinkStart',1,@isnumeric);
            p.addParameter('nSparsePriorPatchStart',0,@isnumeric);
            p.addParameter('specifiedStarts',{},@iscell);
            p.addParameter('sparsePrior',struct([]),@isstruct);

            p.addParameter('maxIter', 1e3, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('bounded', true, @(x)(islogical(x) && numel(x) == 1));
            p.addParameter('ub', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('display', 'iter');
            p.addParameter('gpu', false, @(x)(islogical(x) && numel(x) == 1));

            parse(p, varargin{:});

            % Keep track of all runs
            multistartStruct.runIndex = 0;
            multistartStruct.initImages = {};
            multistartStruct.initTypes = {};
            multistartStruct.reconImages = {};
            multistartStruct.initLosses = []; multistartStruct.initLogPriors = []; multistartStruct.initLogLikelihoods = [];
            multistartStruct.initPreds = [];
            multistartStruct.reconLosses = []; multistartStruct.reconLogPriors = []; multistartStruct.reconLogLikelihoods = [];
            multistartStruct.reconPreds = [];
            multistartStruct.coneVec = coneVec;

            % Some parameters
            meanLuminanceCdPerM2 = [];

            % Run estimates from specified number of white noise starting
            % points
            for ii = 1:p.Results.nWhiteStart
                multistartStruct.runIndex = multistartStruct.runIndex + 1;
                multistartStruct.initTypes{multistartStruct.runIndex} = 'whiteNoise';
                multistartStruct.initImages{multistartStruct.runIndex} = rand(this.Size);
                
                [multistartStruct.reconImages{multistartStruct.runIndex},multistartStruct.initLosses(multistartStruct.runIndex),multistartStruct.reconLosses(multistartStruct.runIndex)] = this.runEstimate(coneVec, ...
                    'init', multistartStruct.initImages{multistartStruct.runIndex}(:), ...
                    'maxIter',p.Results.maxIter,'bounded',p.Results.bounded,'ub',p.Results.ub, ...
                    'display',p.Results.display,'gpu',p.Results.gpu);

                multistartStruct.initPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.initImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.initImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.initLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.initLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.initLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;

                multistartStruct.reconPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.reconImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.reconImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.reconLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.reconLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.reconLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;
            end

            % Run estimates from specified number of pink noise starting
            % points
            for ii = 1:p.Results.nPinkStart
                multistartStruct.runIndex = multistartStruct.runIndex + 1;
                multistartStruct.initTypes{multistartStruct.runIndex} = 'pinkNoise';
                multistartStruct.initImages{multistartStruct.runIndex} = spectrumSampler(this.Size);

                [multistartStruct.reconImages{multistartStruct.runIndex},multistartStruct.initLosses(multistartStruct.runIndex),multistartStruct.reconLosses(multistartStruct.runIndex)] = this.runEstimate(coneVec, ...
                    'init', multistartStruct.initImages{multistartStruct.runIndex}(:), ...
                    'maxIter',p.Results.maxIter,'bounded',p.Results.bounded,'ub',p.Results.ub, ...
                    'display',p.Results.display,'gpu',p.Results.gpu);

                 multistartStruct.initPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.initImages{multistartStruct.runIndex}(:);
                 [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.initImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.initLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.initLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.initLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;

                multistartStruct.reconPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.reconImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.reconImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.reconLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.reconLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.reconLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;
            end

            % Run estimates from specified number of sparse prior patch starting
            % points
            if (p.Results.nSparsePriorPatchStart > 0)
                if (isempty(p.Results.sparsePrior))
                    error('Need to specify sparse prior if using sparse prior patch starts');
                end
            end
            for ii = 1:p.Results.nSparsePriorPatchStart
                multistartStruct.runIndex = multistartStruct.runIndex + 1;
                multistartStruct.initTypes{multistartStruct.runIndex} = 'sparsePriorPatch';
                multistartStruct.initImages{multistartStruct.runIndex} = sparseSampler(p.Results.sparsePrior,this.Size);

                [multistartStruct.reconImages{multistartStruct.runIndex},multistartStruct.initLosses(multistartStruct.runIndex),multistartStruct.reconLosses(multistartStruct.runIndex)] = this.runEstimate(coneVec, ...
                    'init', multistartStruct.initImages{multistartStruct.runIndex}(:), ...
                    'maxIter',p.Results.maxIter,'bounded',p.Results.bounded,'ub',p.Results.ub, ...
                    'display',p.Results.display,'gpu',p.Results.gpu);

                multistartStruct.initPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.initImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.initImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.initLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.initLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.initLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;

                multistartStruct.reconPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.reconImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.reconImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.reconLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.reconLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.reconLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;
            end

            % Run estimates from specified starting points
            for ii = 1:length(p.Results.specifiedStarts)
                multistartStruct.runIndex = multistartStruct.runIndex + 1;
                multistartStruct.initTypes{multistartStruct.runIndex} = 'specified';
                multistartStruct.initImages{multistartStruct.runIndex} = reshape(p.Results.specifiedStarts{ii},this.Size);

                [multistartStruct.reconImages{multistartStruct.runIndex},multistartStruct.initLosses(multistartStruct.runIndex),multistartStruct.reconLosses(multistartStruct.runIndex)] = this.runEstimate(coneVec, ...
                    'init', multistartStruct.initImages{multistartStruct.runIndex}(:), ...
                    'maxIter',p.Results.maxIter,'bounded',p.Results.bounded,'ub',p.Results.ub, ...
                    'display',p.Results.display,'gpu',p.Results.gpu);

                multistartStruct.initPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.initImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.initImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.initLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.initLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.initLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;

                multistartStruct.reconPreds(:,multistartStruct.runIndex) = this.Render * multistartStruct.reconImages{multistartStruct.runIndex}(:);
                [tempNegLogPrior,~,tempNegLogLikely] = ...
                     this.evalEstimate(coneVec, multistartStruct.reconImages{multistartStruct.runIndex}(:));
                if (tempNegLogPrior + tempNegLogLikely ~= multistartStruct.reconLosses(multistartStruct.runIndex))
                    error('Cannot reconstruct loss from prior and likelihood');
                end
                multistartStruct.reconLogPriors(multistartStruct.runIndex) = -tempNegLogPrior;
                multistartStruct.reconLogLikelihoods(multistartStruct.runIndex) = -tempNegLogLikely;
            end

            % Check that we ran at least one estimate
            if (multistartStruct.runIndex == 0)
                error('Need to specify at least one starting scheme');
            end

            % Get best reconstruction out of structure
            [reconImageLinear,reconIndex] = this.selectEstimateFromMultistart(multistartStruct);
        end

        % Choose best reconstruction
        function [reconImageLinear,reconIndex] = selectEstimateFromMultistart(this, multistartStruct)

            % Loop through and find best reconstruction based on loss
            minLoss = Inf;
            reconIndex = NaN;
            for ii = 1:length(multistartStruct.initLosses)
                if (multistartStruct.reconLosses(ii) < minLoss)
                    minLoss = multistartStruct.reconLosses(ii);
                    reconIndex = ii;
                    reconImageLinear = multistartStruct.reconImages{ii};
                end
            end
        end
      
        function [reconstruction,initLoss,solnLoss] = runEstimate(this, coneVec, varargin)
            p = inputParser;
            p.addParameter('maxIter', 1e3, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('init', rand([prod(this.Size), 1]));
            p.addParameter('bounded', true, @(x)(islogical(x) && numel(x) == 1));
            p.addParameter('ub', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('display', 'iter');
            p.addParameter('gpu', false, @(x)(islogical(x) && numel(x) == 1));
            parse(p, varargin{:});

            [reconstruction,initLoss,solnLoss] = this.estimate(coneVec, p.Results.maxIter, p.Results.init, ...
                p.Results.bounded, p.Results.ub, p.Results.display, p.Results.gpu);
        end

        function [reconstruction,initLoss,solnLoss] = estimate(this, measure, maxIter, init, bounded, ub, disp, gpu)
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

            if ~exist('disp', 'var')
                disp = 'iter';
            end

            if ~exist('gpu', 'var')
                gpu = false;
            end

            if gpu
                measure = gpuArray(single(measure));
                loss = @(x) this.reconObjectiveGPU(measure, x);
            end

            % Grab loss to initialization.  Compute return value
            % unnormalized.
            this.LossFactor = 1;
            initLoss = loss(init);
            this.LossFactor = this.LossMultiplier/abs(initLoss);

            if bounded
                    options = optimoptions('fmincon', 'Display', disp, 'MaxIterations', maxIter, 'CheckGradients', false, ...
                    'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, ...
                    'HessianApproximation', 'lbfgs', 'MaxFunctionEvaluations', floor(maxIter * 1.25));
                lb = init * 0;
                ub = ones(size(init)) * ub;
                solution = fmincon(loss, init, [], [], [], [], lb, ub, [], options);
            else
                options  = optimset('GradObj', 'on', 'Display', disp, 'MaxIter', maxIter, 'MaxFunctionEvaluations', floor(maxIter * 1.25));
                solution = fminlbfgs(loss, init, options);
            end

            % Compute returned loss unnormalized.
            this.LossFactor = 1;
            solnLoss = loss(solution);
            reconstruction = reshape(solution, this.Size);
        end

        % This grayscale function does not have all the options that the
        % full estimate function now has.
        function reconstruction = estimateGray(this, measure, maxIter)
            graySize = this.Size([1, 2]); grayLen = prod(graySize);
            rgbSize  = this.Size; rgbLen = prod(rgbSize);

            prmptMtx = zeros(rgbLen, grayLen);
            for idx = 1:rgbLen
                idy = mod(idx, grayLen);
                if idy == 0
                    idy = grayLen;
                end
                prmptMtx(idx, idy) = 1;
            end

            function [loss, gradient] = grayLoss(reconObj, measure, prmptMtx, x)
                [loss, gradient] = reconObj.reconObjective(measure, prmptMtx * x);
                gradient = (gradient' * prmptMtx)';
            end
            
            loss = @(x) grayLoss(this, measure, prmptMtx, x);

            if ~exist('maxIter', 'var')
                maxIter = 1e3;
            end

            init = rand([prod(graySize), 1]);
            this.LossFactor = 1;
            this.InitialLoss = loss(init);

            options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', maxIter, 'CheckGradients', false, ...
                'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, ...
                'HessianApproximation', 'lbfgs', 'MaxFunctionEvaluations', floor(maxIter * 1.25));

            lb = init * 0;
            ub = ones(size(init)) * 1.0;
            solution = fmincon(loss, init, [], [], [], [], lb, ub, [], options);

            this.LossFactor = 1;
            reconstruction = reshape(solution, graySize);
        end

        function [nlogPrior, gradPrior, nlogLlhd, gradLlhd] = evalEstimate(this, measure, imageVec)
            [nlogPrior, gradPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogLlhd,  gradLlhd]  = this.likelihood(measure, imageVec);
            if (~isreal(nlogLlhd))
                error('Log likelihood is not real. That won''t do. Probably there are negative image values.')
            end
        end
    end
end

