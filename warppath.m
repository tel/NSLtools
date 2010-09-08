function [tt, ALIGN] = warppath(ss)
% WARPPATH 1-D warping interpolation (128-point)
%	[tt, ALIGN] = warppath(ss);
%	ss: real positive, M-by-2 matrix
%	tt = t21 + i*t12: warpped axes
%	ALIGN: alignment path
%	Note: This function is not efficient when # > 256

% dimension
if  nargin < 2, lambda = .5; end;
M = size(ss, 1);
L = length(lambda);
sm = zeros(M, L);

% gain
g = max(ss);
for k = 1:2,
	ss(:, k) = ss(:, k) / g(k);
end;

% axis warping
t11 = (1:M)';
t22 = t11;

% distance matrix
d = log(exp(ss(:, 1))*exp(-ss(:, 2).')).^2;

% initiate cumulative distance matrix
DD = d + Inf;
DD(1, 1:2) = cumsum(d(1, 1:2));
DD(1:2, 1) = cumsum(d(1:2, 1));

% initiate link matrix
LL = zeros(M);
LL(1, 1:2) = 0 + j*(0:1);
LL(1:2, 1) = (0:1)';

% cumulation
idx = 0;
for mi = 2:M,
	M1 = max([2, floor(mi/2), 2*mi-M-1]);
	M2 = min([2*mi, ceil((mi+M)/2), M]); 
	for mj = M1:M2,
		if LL(mi-1, mj) == 1,
			[Dij, idx] = min([Inf, DD(mi, mj-1), DD(mi-1, mj-1)]);
		elseif LL(mi, mj-1) == j,
			[Dij, idx] = min([DD(mi-1, mj), Inf, DD(mi-1, mj-1)]);
		else,
			[Dij, idx] = min([DD(mi-1, mj), DD(mi, mj-1), DD(mi-1, mj-1)]);
		end;
		DD(mi, mj) = d(mi, mj) + Dij;
		LL(mi, mj) = rem(idx, 2) + j*floor(idx/2);
	end;
end;

% back trace
ALIGN = zeros(2*M, 1);
ALIGN(1) = M + j*M;
m = 1; mi = M; mj = M;
while mi+mj > 2,
	ALIGN(m+1) = ALIGN(m) - LL(mi, mj);
	mi = real(ALIGN(m+1));
	mj = imag(ALIGN(m+1));
	m = m + 1;
end;
ALIGN = flipud(ALIGN(1:m));

% t2 <-> t1 mapping
t21 = zeros(M, 1);
t12 = zeros(M, 1);
for k = 1:m,
	kr = real(ALIGN(k));
	ki = imag(ALIGN(k));
	t12(kr) = t12(kr) + ki + j;
	t21(ki) = t21(ki) + kr + j;
end;
t21 = real(t21)./imag(t21);
t12 = real(t12)./imag(t12);
tt = [t21+i*t12];

