function image_q(c, c_max, mn)
% IMAGE_Q B&W complex image using quivers
%	image_q(c, c_max, mn)

if nargin < 3, mn = round(size(c)); end;

% gray image
[m, n] = size(c);
CMAP = 1-gray;
C = size(CMAP, 1);
colormap(CMAP);
c_mag = abs(c);
if nargin < 2, c_max = max(c_mag(:)); end;
if isempty(c_max), c_max = max(c_mag(:)); end;

image(c_mag/c_max*C); axis xy;

% sampling
if sum(abs(mn-size(c))),
	mi = linspace(1, m, mn(1));
	ni = linspace(1, n, mn(2));
	ci = interp2(1:n, (1:m)', c, ni, mi') / c_max;
else,
	ci = c;
    mi = 1:m;
    ni = 1:n;
end;
ci = ci * i; 	% rotated by 90 deg
hold on;
quiver(ni, mi, real(ci), imag(ci), .5, 'w-'); 
hold off
