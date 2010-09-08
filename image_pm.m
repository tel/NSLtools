function image_pm(c, c_max, mn, fontsi)
psym = '+'; msym = '-';
if nargin < 2, mn = round(size(c)); end;
if nargin < 3, fontsi = 10; end;

% gray image
[m, n] = size(c);
%CMAP = 1-gray;
CMAP = colormap;
C = size(CMAP, 1);
colormap(CMAP);
c_mag = abs(c);
if nargin < 2, c_max = max(c_mag(:)); end;
if isempty(c_max), c_max = max(c_mag(:)); end;
image(c_mag/c_max*C); axis xy;
%image(c/c_max*C/2+32); axis xy;

% sign symbol
THRESHOLD = .05 * c_max;

if sum(abs(mn-size(c))),
	mi = linspace(1, m, mn(1));
	ni = linspace(1, n, mn(2));
	ci = interp2(1:n, 1:m, c, ni', mi) / c_max;
else,
	ci = c;
    mi = 1:m;
    ni = 1:n;
end;

% interleaving mask
mask_c = (-1).^(1:mn(1));
mask_r = (-1).^(1:mn(2));
mask_m = max(mask_c' * mask_r, 0);
ci = ci .* mask_m;

if length(psym),
    [mdx, ndx] = find(ci>THRESHOLD);
    hp = text(ni(ndx), mi(mdx), psym);
else,
	hp = [];
end;

[mdx, ndx] = find(ci<-THRESHOLD);
hm = text(ni(ndx), mi(mdx), msym);

set([hp; hm], 'fontsi', fontsi, 'co', 'w', 'ho', 'ce', 'clip', 'on');
