function yh = cor2aud2(z, sv, NORM, SRF, M0, FULLOUT)
% COR2AUDS complex inverse wavelet transform.
%	yh = cor2auds(z, sv);
%	yh = cor2auds(z, sv, NORM, SRF, M0, FULLOUT);
%	z: static cortical representation, M-by-K, M = 128.
%	sv: scale vector, K-by-1
%	NORM: (optional) 0=flat, 1=full, .5=partial normalization
%	SRF: (optional) sample ripple frequency (ripple resolution)
%	M0: (optional) original length 
%	FULLOUT: (optional) overlapped output
%
%	COR2AUDS reconstructs auditory spectrum YH from a complex
%	cortical response Z with respect to scale vector SV. 
%	CFIL specify the cortical filters. Normaliztion can be selected
%	out of three style in which FLAT normalization is default. The
%	default ripple resolution is 24 ch / octave ripple.
%	See also: AUD2CORS, COR2AUD

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 20-Aug-97, add non-truncation, M0, FULLOUT
% v1.02: 05-Jan-98, bugs: nargin < 6, dM.
% v1.03: 06-Jan-98, made the 1st (last) filter pure lowpass (highpass)
% v1.04: 12-Apr-98, remove SB option in gen_corf

% target cortical representaion
[M, K] = size(z);
if length(sv) ~= K, error('Size not matched'); end;

% normalization
if nargin < 3, NORM = 0; end;
% sample ripple frequency
if nargin < 4, SRF = 24; end;
if nargin < 5, M0 = M; end;
if nargin < 6, FULLOUT = 0; end;

Z_cum = 0;
HH = 0;
dM = floor((M-M0)/2);
M1 = 2^nextpow2(M0);
M2 = M1 * 2;

if dM,
	z([1:M0+dM (1:dM)+M2-dM], :) = ...
		z([(1:M0+dM)+dM 1:dM], :);
end;

% cumulation
for k = 1:K,
	H	= gen_corf(sv(k), M1, SRF, [k K]);
	R1	= fft(z(:, k), M2);
	Z_cum = Z_cum + R1(1:M1) .* H;
	HH	= HH + H.^2;
end;

% normalization
if nargin < 3,
	NORM = 0;
end;

%HH = HH * NORM + max(HH) * (1 - NORM);
%Z_cum(2:M0) = Z_cum(2:M0) ./ HH(2:M0);

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

if FULLOUT,
	yh = yh([(1:dM)+M2-dM 1:M0+dM]) *2;
else,
	yh = yh(1:M0) *2;
end;
