function h = image_c(c, c_max, mn, fonts, xyaxis)
% IMAGE_C B&W complex image
%	image_c(c, c_max, mn, fontsi, [x0 x y0 y1])

if nargin < 3, mn = round(size(c)); end;
if nargin < 4, fonts = 16; end;
if nargin < 5, xyaxis = [1 mn(2) 1 mn(1)]; end;

UPARROW = setstr(173);
%UPARROW = setstr(221);

% gray image: 3/4 of [w --> k]
[m, n] = size(c);
CMAP = 1 - (linspace(0, 1, 16)'*[1 1 1]*.6);
colormap(CMAP);	% 3/4 of colormap
c_mag = abs(c);
if nargin < 2, c_max = max(c_mag(:)); end;
if isempty(c_max), c_max = max(c_mag(:)); end;

if nargin < 5,
	h = image(c_mag/c_max*32);
else,
	h = image(xyaxis(1:2), xyaxis(3:4), c_mag/c_max*32);
end;
axis xy;

% sign symbol
THRESHOLD = .01;

% sampling
if sum(abs(mn-size(c))),
	mi = linspace(1, m, mn(1));
	ni = linspace(1, n, mn(2));
	ci = interp2(1:n, (1:m)', c, ni, mi');
else,
	ci = c;
    mi = 1:m;
    ni = 1:n;
end;

if nargin >= 5,
	mi = interp1([1 m], xyaxis(3:4), mi);
	ni = interp1([1 n], xyaxis(1:2), ni);
end;

% interleaving mask
mask_c = (-1).^(1:mn(1));
mask_r = (-1).^(1:mn(2));
mask_m = max(mask_c' * mask_r, 0);
ci = ci/c_max .* mask_m;

[mdx, ndx] = find(abs(ci)>THRESHOLD);
for k = 1:length(ndx),
	R1 = abs(ci(mdx(k), ndx(k)));
	if R1 > 1/4, fonts_k = fonts;
	elseif R1 < 1/16, fonts_k = 2;
	else, fonts_k = round((R1-1/16)/(3/16)*(fonts-2))+2; end;
	text('pos', [ni(ndx(k)), mi(mdx(k))], 'str', UPARROW, ...
		'fonts', fonts_k, 'fontn', 'symbol', 'co', 'w', 'ho', 'ce', ...
		'rot', angle(ci(mdx(k), ndx(k)))/pi*180, 'clip', 'on');
end;
