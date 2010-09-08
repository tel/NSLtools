function yh = cor2auds(z, sv, NORM, SRF, M0, FULLOUT, BP)
% COR2AUDS complex inverse wavelet transform.
%	yh = cor2auds(z, sv);
%	yh = cor2auds(z, sv, NORM, SRF, M0, FULLOUT, BP);
%	z: static cortical representation, M-by-K, M = 128.
%	sv: scale vector, K-by-1
%	NORM: (optional) 0=flat, 1=full, .5=partial, (default=.9);
%	SRF: (optional) sample ripple frequency (ripple resolution)
%	M0: (optional) original length 
%	FULLOUT: (optional) overlapped output
%   BP: bandpass filters indicator, default : 0
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

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, bug: nargin < 3 (normalization)
% v1.20: 17-Jan-02, add BP option
% v1.30: 14-Oct-04, add DC normalization for perfect reconstruction

% target cortical representaion
[M, K] = size(z);
if length(sv) ~= K, error('Size not matched'); end;

% normalization
if nargin < 3, NORM = .99; end;
% sample ripple frequency
if nargin < 4, SRF = 24; end;
if nargin < 5, M0 = M; end;
if nargin < 6, FULLOUT = 0; end;
if nargin < 7, BP = 0; end;

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
	H	= gen_corf(sv(k), M1, SRF, [k+BP K+BP*2]);
	R1	= fft(z(:, k), M2);
	Z_cum = Z_cum + R1(1:M1) .* H;
	HH	= HH + H.^2;
end;

% normalization
%if nargin < 3,
%	NORM = 0;
%end;

HH = HH * NORM + max(HH) * (1 - NORM);
HH(1) = 2*HH(1);	 % normalization for DC
Z_cum(1:M0) = Z_cum(1:M0) ./ HH(1:M0);

%plot(abs(Z_cum));
yh = ifft(Z_cum, M2);

if FULLOUT,
	yh = yh([(1:dM)+M2-dM 1:M0+dM]) *2;
else,
	yh = yh(1:M0) *2;
end;
