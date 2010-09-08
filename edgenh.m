function zh = cor_enhs(z, p)
% COR_ENHS edge enhancement algorithm rev. 2.0
%	zh = cor_enhs(z, p);
%	z: an M-by-K matrix;
%	p: enhance mode (multiplied by the negative curvature),
%	p = 1, column-wise enhancement;
%	p = 2, row-wise enhancement (default)
%	p = 3, both. M, N must be greater than 3. 
%	p < 0, amplitude is substituded by its negative curvature.

% Copyright (C) 1995 by P. RU, NSL, ISR, UMD. All rights reserved.
% v1.00: 01-Jan-1995
% v1.01: 03-Jun-1998, abb p < 0 option
[m, n] = size(z);
z0 = z; z = abs(z);
if nargin < 2, p = 2; end;
CURV = (p < 0); p = abs(p);

% column-wise enhancement
if p > 1,	% 2 or 3
	z1 = 0 * z + 1;
else,
	z1 = 0 * z;
	zd = diff(z, 2);
	z1(2:(m-1),:) = zd .* sign(sign(zd)-1);
end;

% row-wise enhancement
if rem(p, 2),	% 1 or 3
	z2 = 0 * z + 1;
else,
	z2 = 0 * z;
	zd = diff(z', 2);
	z2(:,2:(n-1)) = (zd .* sign(sign(zd)-1))';
end;

% total enhancement
zh = (z1 .* z2) .^ (1);
if CURV, zh = zh .* exp(i*angle(z0)); 
else, zh = zh .* z0; end;
