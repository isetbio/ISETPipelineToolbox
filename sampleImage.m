function sample = sampleImage(inputImage, size)
% Sample a small patch of particular size (n by n) from the input image.

nRow = size(inputImage, 1);
nCol = size(inputImage, 2);

idxRow = randi(nRow - n + 1);
idxCol = randi(nCol - n + 1);

sample = inputImage(idxRow:idxRow+n-1, idxCol:idxCol+n-1, :);

end

