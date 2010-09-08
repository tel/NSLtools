function compleximage(varargin)
%
% COMPLEXIMAGE plots a 2 dimensional complex matrix
%
% 	compleximage(matrix)
% 	
% 	compleximage(xaxis,yaxis,matrix)
%
%  The image is plotted with saturation proportional to squared amplitude
%  (where totally unsaturated is gray and totally saturated vivid)
%  and hue proportional to phase (where hue is angle on the RGB color
%  wheel and is periodic). Intensity is the maximum value (1.0 in 
%  Matlab) minus a small constant times the squared amplitude, to
%  provide a little extra contrast.

if nargin == 3
	xaxis  = varargin{1};
	yaxis  = varargin{2};
	matrix = varargin{3};
else
	matrix = varargin{1};
	xaxis = 1:size(matrix,2);
	yaxis = 1:size(matrix,1);
end

xsize = length(xaxis);
ysize = length(yaxis);
rsmat = reshape(matrix,xsize*ysize,1);
rsmat2 = real(rsmat.*conj(rsmat));
rsmat2max = max(max(rsmat2));
rsrgb = hsv2rgb([(angle(rsmat)+pi)/(2*pi),...
			rsmat2/rsmat2max,...
			1-0.2*rsmat2/rsmat2max]);
rgbplot = reshape(rsrgb,ysize,xsize,3);
image(xaxis,yaxis,rgbplot)

