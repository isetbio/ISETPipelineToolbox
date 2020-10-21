function basisImg = visualizeBasis(basisSet, imgDim, basisSize, bw)
% Visualization of grayscale basis function
%
% Syntax:
%   basisImg = visBasisBW(basisSet, imgDim, basisSize)
%
% Description:
%   Visualization of a basis set, given the dimensionality of basis and original images.
%
% Inputs:
%   basisSet   - A matrix of basis functions of size [(imgDim * imgDim), basisSize]
%   imgDim     - The linear dimensionality of original (square) image
%   basisSize  - Number of basis function
%
% Outputs:
%   basisImg   - A single large image consists of each basis function
%
% Examples:
%   imshow(visBasisBW(basisSet, imgDim, basisSize), 'InitialMagnification', 300)

if(bw)
    dx = imgDim;
    dy = imgDim;
    
    basisSet = reshape(basisSet, [dx, dy, basisSize]);
    
    % M by M large "image"
    allDim   = floor(sqrt(basisSize));
    basisImg = zeros(allDim * imgDim, allDim * imgDim);
    
    for i = 1:allDim
        for j = 1:allDim
            % Select each Basis Image
            idx   = (i - 1) * allDim + j;
            basis = basisSet(:, :, idx);
            basis = basis(:);
            
            % Normalization
            basis = (basis - min(basis))/(max(basis) - min(basis));
            
            % Add to Display Image
            basisImg( (i-1) * dx + 1:i * dx, (j-1) * dy + 1:j * dy) = reshape(basis, dx, dy);
        end
    end
else    
    dx = imgDim;
    dy = imgDim;
    
    basisSet = reshape(basisSet, [dx, dy, 3, basisSize]);
    
    % M by M large "image"
    allDim   = floor(sqrt(basisSize));
    basisImg = zeros(allDim * imgDim, allDim * imgDim, 3);
    
    for i = 1:allDim
        for j = 1:allDim
            % Select each Basis Image
            idx   = (i - 1) * allDim + j;
            basis = basisSet(:, :, :, idx);
            basis = basis(:);
            
            % Normalization
            basis = (basis - min(basis))/(max(basis) - min(basis));
            
            % Add to Display Image
            basisImg( (i-1) * dx + 1:i * dx, (j-1) * dy + 1:j * dy, :) = reshape(basis, dx, dy, 3);
        end
    end
end

figure;
imshow(basisImg, 'InitialMagnification', 100);

end