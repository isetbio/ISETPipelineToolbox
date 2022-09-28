function sample = sparseSampler(prior, imSize)
% Sample a large image with patches generated from a sparse image prior
% 
% Syntax:
%   sample = sparseSampler(prior, imSize)
%
% Description
%   Creates a block structured image where each block is drawn from an
%   Laplacian sparse image prior.
%
%   Negative pixel valus are truncated to zero.
%
% Inputs:
%   prior    - Learned sparse prior structure with basis function and mean
%   imSize   - Size of the large image being sampled
%
% Outputs:
%   sample   - Sampled image

% History:
%   09/28/22 dhb  Truncate to non-negative.

stride = sqrt(size(prior.regBasis, 1) / 3);
xStep = ceil(imSize(1) / stride);
yStep = ceil(imSize(2) / stride);

fullImg = zeros(xStep * stride, yStep * stride, 3);

for x = 1 : xStep
    for y = 1 : yStep
        fullImg((x - 1) * stride + 1 : x * stride, ...
                (y - 1) * stride + 1: y * stride, :) ...
                = returnSample(prior);
    end
end

sample = fullImg(1:imSize(1), 1:1:imSize(2), :);
sample(sample < 0) = 0;

end


% Helper function:
% Sample a patch from the sparse image prior
% using an exponential distribution
function sample = returnSample(prior)

stride = sqrt(size(prior.regBasis, 1) / 3);
imSize = [stride, stride, 3];

% sample from an exponential distribution 
rndMu = exprnd(0.3980 * ones(size(prior.mu')));

% sample random binary sign
rndSign = rand(size(rndMu));
rndSign(rndSign > 0.5) = 1; rndSign(rndSign < 0.5) = -1;

% multiply the coefficients with the basis function
sample = prior.mu' + prior.regBasis * (rndMu .* rndSign);

% reshape image
sample = reshape(sample, imSize);
end
