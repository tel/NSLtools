function zij = cor2sing(z)
% COR2SING pit detecting algorithm by phase imparity rev. 2.0
%	zij = cor2sing(z), z is an M-by-N matrix, 
%	Copyright (C) 1995 by P. RU, NSL, ISR, UMD. All rights reserved.

[M, N] = size(z);

% curvature detecting
z = abs(z);
dz_ac = sign(diff(z));
dz_ac = diff(dz_ac); 
dz_ac = max(sign(dz_ac - 1.5), 0);
dz_ac = [zeros(1, N); dz_ac; zeros(1, N)];

dz_ar = sign(diff(z'));
dz_ar = diff(dz_ar);
dz_ar = max(sign(dz_ar - 1.5)', 0);
dz_ar = [zeros(M, 1) dz_ar zeros(M, 1)];

% 
dz = dz_ar.*dz_ac;
[zi, zj, zz] = find(dz);
zij = [zi, zj];

% elimination
zrng = 2;
for k = length(zz):-1:1,
	zref_i = max(1, zi(k)-zrng):min(zi(k)+zrng, M);
	zref_j = max(1, zj(k)-zrng):min(zj(k)+zrng, N);
	zref =  min(min(z(zref_i, zref_j)));
	if z(zi(k), zj(k)) > zref,
		zij(k, :) = [];
	end;
end;

if isempty(zij), zij = [0, 0]; end;
