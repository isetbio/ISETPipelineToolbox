function samples = lgvSampler(prior, nSample, imSize, varargin)
% Implementation of a vanilla version of Langevin sampler for the convolutional sparse prior 

% parameters for the algorithm
p = inputParser;
p.addParameter('burnIn', 250);
p.addParameter('nStep', 250);
p.addParameter('tau', 1e-5);
p.addParameter('gamma', 2.5e-4);
p.addParameter('stride',4)
parse(p, varargin{:});

nStep = p.Results.nStep;
tau = p.Results.tau;
gamma = p.Results.gamma;
stride = p.Results.stride;

% create image estimator object
regPara = 1.0;
estm = PoissonSparseEstimator([], inv(prior.regBasis), ...
                        prior.mu', regPara, stride, imSize);


samples = zeros([nSample, imSize]);

% burn-in
seed = rand(prod(imSize), 1);
for idx = 1:p.Results.burnIn
    [~, grad] = estm.prior(reshape(seed, imSize));
    seed = seed + tau * (-grad) + sqrt(2 * gamma) * normrnd(0, 1, prod(imSize), 1);
end

% take samples
totalStep = nStep * nSample;
sampleIdx = 1;

image = seed;
for idx = 1:totalStep
    [~, grad] = estm.prior(reshape(image, imSize));
    image = image + tau * (-grad) + sqrt(2 * gamma) * normrnd(0, 1, prod(imSize), 1);

    if mod(idx, nStep) == 0  
        samples(sampleIdx, :, :, :) = reshape(image, imSize);
        fprintf('Sample %d / %d \n', sampleIdx, nSample);

        sampleIdx = sampleIdx + 1;      
    end
end

end