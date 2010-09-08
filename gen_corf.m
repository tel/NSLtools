function H = gen_corf(fc, L, SRF, KIND)
% GEN_CORF generate (bandpass) cortical filter transfer function
%	h = gen_corf(fc, L, SRF);
%	h = gen_corf(fc, L, SRF, KIND);
%	fc: characteristic frequency
%	L: length of the filter, power of 2 is preferable.
%	SRF: sample rate.
%	KIND: (scalar)
%	      1 = Gabor function; (optional)
%	      2 = Gaussian Function (Negative Second Derivative) (defualt)
%	      (vector) [idx K]
%	      idx = 1, lowpass; 1<idx<K, bandpass; idx = K, highpass.
%
%	GEN_CORF generate (bandpass) cortical filter for various length and
%	sampling rate. The primary purpose is to generate 2, 4, 8, 16, 32 Hz
%	bandpass filter at sample rate 1000, 500, 250, 125 Hz. This can also
%	be used to generate bandpass spatial filter .25, .5, 1, 2, 4 cyc/oct
%	at sample ripple 20 or 24 ch/oct. Note the filter is complex and
%	non-causal.
%	see also: AUD2COR, COR2AUD
 
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 11-Jun-97, added DSB option.
% v1.02: 28-Aug-97, added KIND, PASS

if nargin < 4, KIND = 2; end;
if length(KIND) == 1,
	PASS = [2 3];	% bandpass
else,
	PASS = KIND;
	KIND = 2;
end;

% fourier transform of lateral inhibitory function 

% tonotopic axis
R1	= (0:L-1)'/L*SRF/2/abs(fc);	% length = L + 1 for now

if KIND == 1,	% Gabor function
	C1      = 1/2/.3/.3;
	H       = exp(-C1*(R1-1).^2) + exp(-C1*(R1+1).^2);
else,		% Gaussian Function 
	R1	= R1 .^ 2;		% 
	H	= R1 .* exp(1-R1); 	% single-side filter
	%H	= H .^ .25;
end;
%H = H * (1-exp(-(fc/.5).^0.5));

% passband
if PASS(1) == 1,		%lowpass
	[maxH, maxi] = max(H);
	sumH = sum(H);
	H(1:maxi-1) = ones(maxi-1, 1);
	H = H / sum(H) * sumH;
elseif PASS(1) == PASS(2),	% highpass
	[maxH, maxi] = max(H);
	sumH = sum(H);
	H(maxi+1:L) = ones(L-maxi, 1);
	H = H / sum(H) * sumH;
end;
