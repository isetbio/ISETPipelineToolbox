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

        % use GPU (gpuArray) to compute the large matrix product in the
        % likelihood function calculation
        function [loss, gradient] = reconObjectiveGPU(this, measure, imageVec)
            [nlogllPrior, gradientPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogllLlhd,  gradientLlhd]  = this.likelihood(measure, gpuArray(single(imageVec)));

            loss = nlogllPrior + double(gather(nlogllLlhd));
            gradient = gradientPrior + double(gather(gradientLlhd));
        end

        function [reconstruction,initLoss,solnLoss] = runMultistartEstimate(this, coneVec, varargin)
            p = inputParser;
            p.KeepUnmatched = false;
            p.addParameter('nWhiteStart', 0, @isnumeric);
            p.addParameter('nPinkStart',1,@isnumeric);
            p.addParameter('nPatchStart',1,@isnumeric);
            p.addParameter('specifiedStarts',{},@iscell);

            p.addParameter('maxIter', 1e3, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('bounded', true, @(x)(islogical(x) && numel(x) == 1));
            p.addParameter('ub', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('display', 'iter');
            p.addParameter('gpu', false, @(x)(islogical(x) && numel(x) == 1));

            parse(p, varargin{:});

            % Keep track of all runs
            runIndex = 0;
            initImages = {};

            % Run estimates from specified number of white noise starting
            % points
            for ii = 1:nWhiteStart
                runIndex = runIndex + 1;
                initImages{runIndex} = rand([prod(this.Size));
                [reconstructions{runIndex},initLosses(runIndex),solnLosses(runIndex)] = this.runEstimate(coneVec, ...
                    'init', initImages{runIndex}, ...
                    'maxIter',p.Results.maxIter,'bounded',p.Results.bounded,'ub',p.Results.ub, ...
                    'display',p.Results.display,'gpu',p.Results.gpu);
            end

            % Check that we ran at least one estimate
            if (runIndex == 0)
                error('Need to specify at least one starting scheme');
            end
           
        end


        % meanLuminanceCdPerM2 = [];
        % scaleFactor = (forwardPupilDiamMM/reconPupilDiamMM)^2;
        % [recon1Image,recon1InitLoss,recon1SolnLoss] = estimator.runEstimate(forwardExcitationsToStimulusUse * scaleFactor, ...
        %     'maxIter', 500, 'display', 'iter', 'gpu', false, 'init', 0.5*ones(length(stimulusImageLinear(:)), 1));
        % [recon1Scene, ~, recon1ImageLinear] = sceneFromFile(gammaCorrection(recon1Image, theForwardDisplay), 'rgb', ...
        %     meanLuminanceCdPerM2, forwardConeMosaic.Display);
        % recon1Scene = sceneSet(recon1Scene, 'fov', fieldSizeDegs);
        % [recon1NegLogPrior,~,recon1NegLogLikely] = ...
        %     estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, recon1ImageLinear(:));
        % visualizeScene(recon1Scene, 'displayRadianceMaps', false,'avoidAutomaticRGBscaling', true);
        % saveas(gcf,fullfile(outputDir,'Recon1.jpg'),'jpg');

        % Report back the better
        % if (-(recon1NegLogPrior+recon1NegLogLikely) > -(recon2NegLogPrior+recon2NegLogLikely))
        %     reconWhichStr = 'Recon1 (random start) better\n';
        %     reconNegLogPrior = recon1NegLogPrior;
        %     reconNegLogLikely = recon1NegLogLikely;
        %     reconImage = recon1Image;
        %     reconScene = recon1Scene;
        %     reconImageLinear = recon1ImageLinear;
        % else
        %     reconWhichStr = 'Recon2 (stimulus start) better\n';
        %     reconNegLogPrior = recon2NegLogPrior;
        %     reconNegLogLikely = recon2NegLogLikely;
        %     reconImage = recon2Image;
        %     reconScene = recon2Scene;
        %     reconImageLinear = recon2ImageLinear;
        % end
        % 
        % % Show reconstruction
        % [reconScene, ~, reconImageLinear] = sceneFromFile(gammaCorrection(recon2Image, theForwardDisplay), 'rgb', ...
        %     meanLuminanceCdPerM2, forwardConeMosaic.Display);
        % reconScene = sceneSet(reconScene, 'fov', fieldSizeDegs);
        % visualizeScene(reconScene, 'displayRadianceMaps', false, 'avoidAutomaticRGBscaling', true);
        % saveas(gcf,fullfile(outputDir,'Recon.jpg'),'jpg');
        % 
        % % Compute forward excitations from reconstruction
        % % And compare with stimulus exciations
        % forwardOI = oiCompute(reconScene,forwardOI);
        % if (reconstructfromRenderMatrix)
        %     title('Reconstruction from forward render matrix');
        %     forwardExcitationsToRecon = squeeze(forwardRenderMatrix*reconImageLinear(:));
        % else
        %     title('Reconstruction from forward ISETBio');
        %     forwardExcitationsToRecon = squeeze(forwardConeMosaic.Mosaic.compute(forwardOI, 'opticalImagePositionDegs', 'mosaic-centered'));
        % end
        % figure; clf; hold on;
        % plot(forwardExcitationsToStimulusUse,forwardExcitationsToRecon,'ro','MarkerFaceColor','r','MarkerSize',10);
        % axis('square');
        % maxVal = max([forwardExcitationsToStimulusUse; forwardExcitationsToRecon]);
        % plot([0 maxVal],[0 maxVal],'k');
        % xlim([0 maxVal]); ylim([0 maxVal]);
        % xlabel('Excitations to stimulus');
        % ylabel('Excitations to reconstruction');
        % saveas(gcf,fullfile(outputDir,'StimulusVsReconExcitations.jpg'),'jpg');
        % 
        % %% Evaluate prior and likelihood of stimulus and reconstruction
        % [stimNegLogPrior,~,stimNegLogLikely] = ...
        %     estimator.evalEstimate(forwardExcitationsToStimulusUse * scaleFactor, stimulusImageLinear(:));
        % txtFileName = fullfile(outputDir,'ReconProbInfo.txt');
        % if (exist(txtFileName,'file'))
        %     delete(txtFileName);
        % end
        % fid = fopen(txtFileName,'w');
        % fprintf(fid,'Stimulus: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
        %     -stimNegLogPrior,-stimNegLogLikely,-(stimNegLogPrior+stimNegLogLikely));
        % fprintf(fid,'Recon1: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
        %     -recon1NegLogPrior,-recon1NegLogLikely,-(recon1NegLogPrior+recon1NegLogLikely));
        % fprintf(fid,'Recon1 initial loss %0.6g; recon1 solution loss %0.6g; fractional difference (init less soln; should be pos): %0.6g\n', ...
        %     recon1InitLoss,recon1SolnLoss,(recon1InitLoss-recon1SolnLoss)/abs(recon1InitLoss));
        % fprintf(fid,'Recon2: reg weighted log prior %0.6g; estimate part of log likelihood %0.6g; sum %0.6g\n', ...
        %     -recon2NegLogPrior,-recon2NegLogLikely,-(recon2NegLogPrior+recon2NegLogLikely));
        % fprintf(fid,'Recon2 initial loss %0.6g; recon2 solution loss %0.6g; fractional difference (init less soln; should be pos): %0.6g\n', ...
        %     recon2InitLoss,recon2SolnLoss,(recon2InitLoss-recon2SolnLoss)/abs(recon2InitLoss));
        % fprintf(fid,reconWhichStr);
        % fprintf(fid,'Each of the following should be *higher* for a valid reconstruction\n');
        % if (-stimNegLogPrior > -reconNegLogPrior)
        %     fprintf(fid,'\tReconstruction prior *lower* than stimulus\n');
        % else
        %     fprintf(fid,'\tReconstruction prior *higher* than stimulus\n');
        % end
        % if (-stimNegLogLikely > -reconNegLogLikely)
        %     fprintf(fid,'\tStimulus likelihood *higher* than reconstruction\n');
        % else
        %     fprintf(fid,'\tStimulus likelihood *lower* than reconstruction\n');
        % end
        % if (-(stimNegLogPrior+stimNegLogLikely) > -(reconNegLogPrior+reconNegLogLikely))
        %     fprintf(fid,'\tReconstruction neg objective *lower* than stimulus\n');
        % else
        %     fprintf(fid,'\tReconstruction neg objective *higher* than stimulus\n');
        % end
        % fclose(fid);

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

            initLoss = loss(init);
            solnLoss = loss(solution);
            reconstruction = reshape(solution, this.Size);
        end

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
            options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', maxIter, 'CheckGradients', false, ...
                'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, ...
                'HessianApproximation', 'lbfgs', 'MaxFunctionEvaluations', floor(maxIter * 1.25));

            lb = init * 0;
            ub = ones(size(init)) * 1.0;
            solution = fmincon(loss, init, [], [], [], [], lb, ub, [], options);

            reconstruction = reshape(solution, graySize);
        end

        function [nlogPrior, gradPrior, nlogLlhd, gradLlhd] = evalEstimate(this, measure, imageVec)
            [nlogPrior, gradPrior] = this.prior(reshape(imageVec, this.Size));
            [nlogLlhd,  gradLlhd]  = this.likelihood(measure, imageVec);
        end
    end
end

