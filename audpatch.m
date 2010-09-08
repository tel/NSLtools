function yh = audpatch(y1, y2, sv, m)
% AUDPATCH auditory patch
%	yh = audpatch(y1, y2, sv, m);
%	y1, y2: auditory spectrogram, N-by-M.
%	sv: scale vector
%	This function basically mixes the "timbre" of [y1] and the pitch or melody 
%	of [y2] to generate a new sound.
%	if they are different in length, the shorter one will be interpolated.
%	The temporal energy of [y1] will be considered. 

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, add upper bound for GAP

% interpolation
[N1, M] = size(y1);
[N2, M] = size(y2);
if N1 == N2,
	N = N1;
else,
	N = sqrt(N1*N2);
	y2 = interp1(1:N2, y2, linspace(1, N2, N));
	y1 = interp1(1:N1, y1, linspace(1, N1, N));
end;

% scale mask
K = length(sv);
GAP = 3;
%FILT = 6;
%B = cos((-(FILT-1):(FILT-1))/FILT*pi/2).^2;
if length(m) == 1,
	m1 = [ones(M, m) zeros(M, K-m)];
	m2 = 1-m1;
else,
	% m(1), m(2)
	m = round(linspace(m(1), m(2), K));
	GAP = min(GAP, m(1));
	m1 = zeros(M, K);
	for k = 1:K,
		m1(:, k) = [zeros(m(k)+GAP,1); ones(M-m(k)-GAP,1)];
		m2(:, k) = [ones(m(k)-GAP,1); zeros(M-m(k)+GAP,1)];
	end;
end;

t0 = clock;
for n = 1:N,
	E1 = sum(y1(n, :).^2);
	E2 = sum(y2(n, :).^2);
	y2(n, :) = sqrt(E1/E2)*y2(n, :);
	% y -> z
	z1 = aud2cors(log(y1(n, :)), sv);
	z2 = aud2cors(log(y2(n, :)), sv);

	% manipulation
	zh = m1 .* z1 + m2 .* z2;
		
	% z -> y
	yh(n, :) = cor2auds(zh, sv).'; 
	time_est(n, N, 50, t0);
end;
yh = exp(real(yh));
