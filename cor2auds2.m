function yh = cor2auds(z, sv, NORM, SRF)
% COR2AUDS complex inverse wavelet transform.
%	XH = A1WT_INV(YY, CFIL) is the reconstructed auditory spectrum XH for
%	a full complex cortical response YY. 
%	CFIL specify the cortical filters.
%
%	Because MATLAB only takes two dimensional matrix, X must be a vector.
%	(c) 1993, 1994 Kuansan Wang, all rights reserved.
%       (c) 1995, 1996 Powen Ru, all rights reserved.
%
%	filter bank:   	 cfil = M-by-N matrix
%	desired output:  yy   = M-by-(N+L-1) matrix
%	estimate input:	 xh   = 1-by-L vector

% sample ripple frequency
if nargin < 4, SRF = 24; end;

% target cortical representaion
[M, K] = size(z);
if length(sv) ~= K, error('Size not matched'); end;
M1 = 2^nextpow2(M);
M2 = M1 * 2;

Z_cum = 0;
HH = 0;

% cumulation
for k = 1:K,
	H	= gen_corf(sv(k), M1, SRF);
	R1	= fft(z(:, k), M2);
	Z_cum	= Z_cum + R1(1:M1) .* H;
	HH	= HH + H.^2;
end;

% normalization
if nargin < 3,
	NORM = 0;
end;

if NORM == 1,
	m1 = 35;
	HH(m1:M) = HH(m1:M) * 0 + HH(m1);
	Z_cum(2:M) = Z_cum(2:M) ./ HH(2:M);
elseif NORM == 0,
	Z_cum = Z_cum / max(HH);
else,
	[maxH, idx] = max(HH);
	HH(idx:M) = HH(idx:M) * maxH;
	Z_cum(2:M) = Z_cum(2:M) ./ HH(2:M);
end;

%plot(abs(Z_cum));
yh = ifft(Z_cum, M2);
yh = yh(1:M) *2;
