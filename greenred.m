function b = greenred(m)
%GREENRED   Bilinear green/red color map.
%
%   GREENRED(M) returns an M-by-3 matrix containing the colormap.
%   GREENRED, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(greenred)
%
%
%   See also REDBLUE, BROWNBLUE, YARG, YARGPRINT, GRAYPRINT.
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
top = floor(m/2); 
tt = [1:top]'/top; % i.e. ('half' of m) numbers from 0 to 1
bot = m-top; 
bb = [bot:-1:1]'/bot; % i.e. ('other half' of m) numbers from 1 to 0

% Hue: 0.35*ones = good green
% Saturation: tt.^2 = graded saturation, more (than linear) at end of scale
% Value: (1-0.3*tt) = graded value from 1 (bright) to 0.3 less bright
btop = [0.35*ones(top,1),tt.^2,(1-0.5*tt)];
% 1.0*ones hue = good red
% Hue: 1.0*ones = good red
% Saturation: (1-bb).^2= graded saturation, more (than linear) at end of scale
% Value: ones = graded value from 1 (bright) to 0.1 barely less bright
bbot = [1.0*ones(bot,1),bb.^2,(1-0.2*bb)];
b=hsv2rgb([bbot;btop]);

if 0
	btop = [1.0*ones(top,1),([1:top]'/top).^2,ones(top,1)];   % 1.0 hue good red
	bbot = [0.7*ones(bot,1),(1-[1:bot]'/bot).^2,ones(bot,1)]; % 0.7 hue good blue
	b=hsv2rgb([bbot;btop]);
end
