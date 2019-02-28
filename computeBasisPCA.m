function [regBasis, mu] = computeBasisPCA(imageSet, imageDim)
%COMPUTEBASISPCA

[pcaBasis, ~, pcaVar, ~, ~, mu] = pca(imageSet);

% Visulization
basisSize = imageDim * imageDim * 3;
basisSet = reshape(pcaBasis, [imageDim, imageDim, 3, basisSize]);

allDim   = floor(sqrt(basisSize)); 
basisImg = zeros(allDim * imageDim, allDim * imageDim, 3);

for i = 1:allDim
    for j = 1:allDim
        % Select each Basis Image
        idx   = (i - 1) * allDim + j;
        basis = basisSet(:, :, :, idx);
        basis = basis(:);
        
        % Normalization
        basis = (basis - min(basis))/(max(basis) - min(basis));

        % Add to Display Image 
        basisImg( (i-1) * imageDim + 1:i * imageDim, (j-1) * imageDim + 1:j * imageDim, :) ...
            = reshape(basis, imageDim, imageDim, 3);
    end
end

figure;
imshow(basisImg, 'InitialMagnification', 50);

scaleMatrix = diag(sqrt(pcaVar));
regBasis    = pcaBasis * scaleMatrix;

end

