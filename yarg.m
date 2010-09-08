function g = yarg(m)
%YARG   Linear gray-scale color map, with opposite dirction of 'gray'.
%   YARG(M) returns an M-by-3 matrix containing a gray-scale colormap.
%   YARG, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(yarg)
%
%   See also REDBLUE, BROWNBLUE, YARGPRINT, GRAYPRINT.
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
g = (0:m-1)'/max(m-1,1);
g = 1-[g g g];
