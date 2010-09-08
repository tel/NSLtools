function yh = audskew(y, sv, UP)
% AUDPATCH 
%
if nargin < 3, UP = 1; end;
[N, M] = size(y);
K = length(sv);

t0 = clock;
for n = 1:N,
	z = aud2cors(y(n, :), sv);
	if UP,
		for k = 2:K, z(:, k) = [zeros(k-1, 1); z(1:(M-k+1), k)]; end;
	else,
		for k = 2:K, z(:, k) = [z(k:M, k); zeros(k-1, 1)]; end;
	end;
	yh(n, :) = cor2auds(z, sv).'; 
	time_est(n, N, 4, t0);
end;
