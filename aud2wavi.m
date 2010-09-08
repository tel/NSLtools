function x0 = aud2wavi(v5, paras)
% AUD2WAVI initial guess for fast inverse auditory spectrum
%	x0 = aud2wav1(v5, [L_frm, tc, fac, shft]);
%	v5		: auditory spectrogram (N-by-M)
%	x0		: the guessed acoustic input.
%
%	PARAS	= [L_frm, tc, fac, shft];
%	L_frm	: frame length, typically, 16 ms or 2^[natural #] ms.
%	tc		: time const, typically, 64 ms = 1024 pts for 16 kHz.
%			if tc == 0, the leaky integration turns to short-term average
%	fac		: nonlinear factor (critical level ratio), typically, .01.
%			The less the value, the more the compression
%			fac = 0: y = (x > 0), full compression
%			fac = -1, y = max(x, 0), half-wave rectifier
%			fac = -2, y = x, linear function
%	shft	: shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
%			etc. SF = 16K * 2^[shft].
%
%	AUD2WAVI estimate the initial guess for fast inverse auditory spectrum.
%	See also: WAV2AUD	

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 20-Jan-98
[N, M] = size(v5);
x0 = 0;
SF = 16000 * 2^paras(4);

CF = cochfil(1:128, paras(4))/SF;
T = paras(1)*SF/1000;	% # of points per frame
L = N * T;				% Total # of points

v5 = v5(:, 1:48);

for m = 1:48,
	x = cos(2*pi*(CF(m)*(1:L)))';
	m = ones(T, 1) * v5(:, m)';
	x = m(:) .* x;
	x0 = x0 + x;
end;
x0 = unitseq(x0);
