function x = sg2sig(H, dT, fr, SF, pg)
% SG2SIG spectrogram to signal (direct synthesis)
%       x = sg2sig(H, dT, fr, SF, pg)
%	H:  M-channel frequency response, (N-by-M)
%	dT: time resolution, (ms) 
%	fr: freq. range vector. (Hz)
%	SF: sample freq. (Hz)
%       pg: power gain,
%
%	SG2SIG synthesizes time-varying multi-tone signal directly from the
%	amplitude spectrogram. The starting phase for each channel is assigned
%	randomly. The sound quality may not be fine. This could be used to
%	compute the initial guess for some spectrogram reconstruction
%	algorithm.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 06-Jun-94

[N, M] = size(H);

% modulated by power gain
if nargin > 4,
	for n = 1:N,
		H(n, :) = normr(H(n, :)) .* pg;
	end;
end;

% define frequency, time and phase
omega	= 2*pi*fr(:)/SF;
phi		= rand(M, 1) * 2 * pi;
dT		= round(dT*SF/1000);
NdT		= N*dT;

t	= (1:NdT)';
LAP	= round(2^19/NdT);
t0	= clock;
x	= 0;

for m = 1:M,
	Hi	= ones(dT, 1) * H(:, m)';
	tones	= cos(omega(m)*t + phi(m));
	x	= x + Hi(:) .* tones;
	time_est(m, M, LAP, t0);
end;

% lowpass filter
%[b, a]	= butter(4, .3);
b = [0.0186    0.0743    0.1114    0.0743    0.0186];
a = [1.0000   -1.5704    1.2756   -0.4844    0.0762];
%x	= filter(b, a, x);

% 50 ms onset or less
rt	= min(fix(.05 * SF), length(x)); 
x(1:rt) = x(1:rt) .* sin(pi/2*(1:rt)/rt)';
