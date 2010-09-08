function h = imagesc0(varargin)
%IMAGESC0 Scale data and display as image, centering on zero.
%   IMAGESC0(...) is the same as IMAGESC(...) except the data is scaled
%   to use to center the colormap on zero.

imagetoplot = varargin{end};
maxsc = max([max(max(imagetoplot)) -min(min(imagetoplot))]);
newvararg = varargin;
newvararg{length(varargin)+1} = [-1 1]*maxsc;
hh = imagesc(newvararg{:});

if nargout > 0
    h = hh;
end
