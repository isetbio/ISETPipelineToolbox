function theDisplay = ScaleDisplayPrimaries(theDisplay, scaleVector, varargin)
% Display version where each rgb channel of monitor
% representation has been scaled by the corresponding entry
% of scaleVector. Used to tweak monitor scaling without needing
% to recompute all of the nice render matrices we may already have.
%
% See also ScaleRenderMatrix

forwardPrimaries = displayGet(theDisplay,'spd primaries')*diag(scaleVector);
theDisplay = displaySet(theDisplay,'spd primaries',forwardPrimaries);

end