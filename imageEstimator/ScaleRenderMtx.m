function renderMtx = ScaleRenderMtx(renderMtx, scaleVector, varargin)
% Render matrix version where each rgb channel of monitor
% representation has been scaled by the corresponding entry
% of scaleVector. Used to tweak monitor scaling without needing
% to recompute all of the nice render matrices we may already have.
%
% See also ScaleDisplayPrimaries

[~,nCols] = size(renderMtx);
nPrimaries = length(scaleVector);
nPixels = nCols/nPrimaries;
if (nPixels ~= floor(nPixels))
    error('Logic error recreating number of pixels from render matrix');
end

for jj = 1:nPixels
    for pp = 1:nPrimaries
        renderMtx(:,jj+(pp-1)*nPixels) = renderMtx(:,jj+(pp-1)*nPixels)*scaleVector(pp);
    end
end

end