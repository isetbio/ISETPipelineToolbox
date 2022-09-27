function sample = spectrumSampler(imSize)
% Sample a large image with 1/f specturm
% 
% Syntax:
%   sample = spectrumSampler(imSize)
%
% Inputs:
%   imSize   - Size of the image being sampled
%
% Outputs:
%   sample   - Sampled image

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

% Transform back to color space
sample = reshape(sample, imSize(1) * imSize(2), imSize(3));
sample = (pca_mtx \ sample')';

% Normalize the image
sample = sample - min(sample(:));
sample = sample / max(sample(:));

sample = reshape(sample, imSize);

end

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