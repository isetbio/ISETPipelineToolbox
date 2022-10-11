function sample = spectrumSampler(imSize, gain)
% Sample a large image with 1/f specturm
% 
% Syntax:
%   sample = spectrumSampler(imSize)
%
% Description:
%   Simple 1/f noise color image sampler.  Draws pink noise in three
%   roughly decorrelated color planes and transforms back to correlated
%   image space.
%
%   Returned pixel values are forced non-negative by truncation to zero.
%
% Inputs:
%   imSize   - Size of the image being sampled
%   gain     - (optional) Pass the pink noise image through a nonlinear
%              function that produces edge structures
%
% Outputs:
%   sample   - Sampled image

% History:
%   09/28/22 dhb  Truncate to non-negative.

% default gain is empty
if ~exist('gain', 'var')
    gain = [];
end

% Define the transformation from pixel to uncorrected color space
pca_mtx = [0.57, 0.5825, 0.5714;
          -0.68, -0.0378, 0.7295; 
          -0.4465, 0.8119, -0.376];

std_cmp = sqrt([0.2169, 0.0175, 0.0044]);

% Sample 1/f noise for each plane
sample = zeros(imSize);
for idx = 1:3
    sample(:, :, idx) = returnSample(imSize(1:2)) * std_cmp(idx);
end

% Pass through a nonlinear function
if ~isempty(gain)    
    phi = erf(gain * sample / sqrt(2));
    z = 2 / pi * asin(gain ^ 2 / (1 + gain ^ 2));
    sample = phi / sqrt(z);
end

% Transform back to color space
sample = reshape(sample, imSize(1) * imSize(2), imSize(3));
sample = (pca_mtx \ sample')';

% Normalize the image
sample = sample - min(sample(:));
sample = sample / max(sample(:));

sample = reshape(sample, imSize);

sample(sample < 0) = 0;
end

% Generate an image with 1/f specturm with random phase
function image = returnSample(imSize)

imCenter = floor(imSize / 2);
amplitude = zeros(imSize);
phase = exp(1i * (rand(imSize) * 2 * pi));

for x = 1:imSize(1)
    for y = 1:imSize(2)
        amplitude(x, y) = ...
            1 / sqrt((x - imCenter(1)) ^ 2 + ...
                     (y - imCenter(2)) ^ 2 + 1);
    end
end

specturm = 1e3 * amplitude .* phase;
image = real(ifft2(ifftshift(specturm)));

end