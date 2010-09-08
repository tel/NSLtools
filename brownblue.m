function b = brownblue(m)
%BROWNBLUE   Bilinear brown/blue color map.
%
%   BROWNBLUE(M) returns an M-by-3 matrix containing the colormap.
%   BROWNBLUE, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(brownblue)
%
%   See also REDBLUE, YARG, YARGPRINT, GRAYPRINT.
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
top = floor(m/2);
bot = m-top;

% btop = interp1([0;1],[.36,.24,.14;1,1,1],([1:top]'/top).^0.5);
btop = interp1([0;1],[.5,.24,.14;1,1,1],([1:top]'/top).^0.5);
bbot = interp1([0;1],[0,0,1;1,1,1],([1:bot]'/bot).^0.4);
b = [bbot; flipud(btop)];
