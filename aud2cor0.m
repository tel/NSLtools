function z = aud2cor0(y, sv, SRF, FULL, BP)
% AUD2CORS static cortical representation.
%	z = aud2cors(y, sv);
%	z = aud2cors(y, sv, SRF, FULL, BP);
%	y	: auditory spectrum (M-by-1)
%	sv	: char. ripple freq's (K-by-1), e.g., sv = 2.^(-2:.5:3)
%	SRF	: (optional) sample ripple frequency, default = 24 ch/oct
%	z	: cortical representation (M-by-K)
%	FULL: non-truncation factor
%	BP	: all bandpass filter bank
%
%	AUD2CORS computes the K-channel cortical representation for an 
%	auditory spectrum Y. The ANGLE(Z) represents the symmetry the unit
%	with the maximum response ABS(Z). For display purpose, the complex
%	matrix Z should be encoded to the indice of a specific colomap by
%	an utility program CPLX_COL.
%	See also: WAV2AUD, AUD2CORS, CPLX_COL, COR_DIST, COR_MAPS

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 19-Aug-97, new: FULL, dM.
% v1.02: 06-Jan-98, made the 1st (last) filter pure lowpass (highpass)
% v1.03: 12-Apr-98, remove SB option in gen_corf
% dimension and zero-padding
K = length(sv);
M = length(y);
M1 = 2^nextpow2(M);
M2 = M1 * 2;
Mp = round((M2-M)/2);

if nargin < 3, SRF = 24; end;	% 24 channel per octave
if nargin < 4, FULL = 0; end;	% with truncation
if nargin < 5, BP = 0; end;

% fourier transform
ypad = y(M) + (1:(M2-M))/(M2-M+1) * (y(1) - y(M)); 
%ypad = [y(M)*ones(Mp, 1); y(1)*ones(M1-Mp, 1)];

y = [y(:); ypad(:)];
Y = fft(y(:));
Y = Y(1:M1);

% spatial filtering
dM = floor(M/2*FULL);
if FULL,
	mdx = [(1:dM)+M2-dM 1:M+dM];
else,
	mdx = 1:M;
end;

z = zeros(M+2*dM, K);
for k = 1:K,
	H	= gen_corf(sv(k), M1, SRF, [k+BP K+BP*2]);
	R1	= ifft(H.* Y, M2);
	z(:, k) = R1(mdx);
end;

% find the energy center --> dilation center
a = abs(z);
p = angle(z);
a_mean = mean(a');
x0 = a_mean*(1:M)'/sum(a_mean);

% warpping
for k = 1:K,
	x_axis = ((1:M)-x0)/(sv(k)) + x0;
	if sv(k) > 1, 		% expand it == interpolation
		a(:, k) = interp1(1:M, a(:, k), x_axis)';
		p(:, k) = interp1(1:M, p(:, k), x_axis)';
	elseif sv(k) < 1,	% shrink it == extrapolation
		x0_axis = [x_axis(1) 1:M, x_axis(M)];
		p1 = (p(1)-p(2))*(1-x_axis(1))+p(1);
		pM = (p(M)-p(M-1))*(x_axis(M)-M)+p(M);
		a(:, k) = interp1(x0_axis, [0; a(:, k); 0], x_axis)';
		p(:, k) = interp1(x0_axis, [p1; p(:, k); pM], x_axis)';
	end;
	%a(:, k) = a(:, k) / sv(k);
end;

z = a .* exp(i*p);
