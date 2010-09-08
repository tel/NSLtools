function H = gen_cort(fc, L, STF, PASS)
% GEN_CORT generate (bandpass) cortical temporal filter transfer function
%	h = gen_cort(fc, L, STF);
%	h = gen_cort(fc, L, STF, PASS);
%	fc: characteristic frequency
%	L: length of the filter, power of 2 is preferable.
%	STF: sample rate.
%	PASS: (vector) [idx K]
%	      idx = 1, lowpass; 1<idx<K, bandpass; idx = K, highpass.
%
%	GEN_CORT generate (bandpass) cortical temporal filter for various
%	length and sampling rate. The primary purpose is to generate 2, 4,
%	8, 16, 32 Hz bandpass filter at sample rate ranges from, roughly
%	speaking, 65 -- 1000 Hz. Note the filter is complex and non-causal.
%	see also: AUD2COR, COR2AUD, MINPHASE
 
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 11-Jun-97, added DSB option.
% v1.02: 28-Aug-97, added KIND, PASS

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, change impulse response of filter 

if nargin < 4, PASS = [2 3]; end;


% tonotopic axis
%R1	= (0:L)'/L*STF/2;	% length = L + 1 for now
%R2	= 15; 
% Gamma distribution function 
% H(f) = f^alpha exp(-beta f)
%H = (R1/fc).^(R2) .* exp(-R2*(R1/fc));
t = (0:L-1)'/STF * fc;
%h = cos(2*pi*t) .* sqrt(t) .* exp(-2*t) * fc;
%h = cos(2*pi*t) .*t.^3 .* exp(-4*t) * fc;
h = sin(2*pi*t) .*t.^2.* exp(-3.5*t) * fc;

%h = diff(h); h = [h(1)/2; h];
h = h-mean(h);
H0 = fft(h, 2*L);
A = angle(H0(1:L));
H = abs(H0(1:L));
[maxH, maxi] = max(H);
H = H / max(H);

% passband
if PASS(1) == 1,		%lowpass
	%sumH = sum(H);
	H(1:maxi-1) = ones(maxi-1, 1);
	%H = H / sum(H) * sumH;
elseif PASS(1) == PASS(2),	% highpass
	%sumH = sum(H);
	H(maxi+1:L) = ones(L-maxi, 1);
	%H = H / sum(H) * sumH;
end;

H = H .* exp(i*A);

