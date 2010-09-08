function b = yargprint(m,cmp,zer)
%YARGPRINT   Nonlinear grayscale colormap.
%
%   YARGPRINT(M) returns an M-by-3 matrix containing the colormap.
%   YARGPRINT, by itself, is the same length as the current colormap.
%   YARGPRINT(M, COMPRESSION, ZEROSHIFT) allows different values for
%     COMPRESSION of the scale near the black/white extremes (default 0.7, 
%     allowed values > 0 ) and ZEROSHIFT of the gray value which
%     represents the center of the black/white continuum (default 0.4,
%     allowed values > 0).
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(yargprint)
%

if nargin < 1, m = size(get(gcf,'colormap'),1); end
if ~exist('cmp','var'); cmp = .70; end
if ~exist('zer','var'); zer = .40; end

gf=gray(m)-.5;
b = ((1+sin(pi*sign(-gf).*((2.*(abs(gf))).^cmp)/2))/2).^zer;

% The sin() changes a line into the sigmoidal part of a sinuoisoid 
% (e.g. from -pi/2 to pi/2). The cmp parameter distorts the line before
% making the sigmod, effectively controlling the slope at the center
% of the sigmoid. The zer parameter then distorts the sigmoid, allowing the
% center of the sigmoid to have a different value than usual (since that value
% usually looks darker in print than halfway between black and white).
