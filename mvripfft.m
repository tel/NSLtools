function [s, ph_c, fdx] = mvripfft(para, cond, ph_c);
% MVRIPFFT generate a single moving ripple via FFT
%	s = mvripfft(para);
%	[s, ph_c, fdx] = mvripfft(para, cond, ph_c);
%	para = [Am, Rt, Om, Ph];
%		Am: modulation depth, 0 < Am < 1, DEFAULT = .9;
%		Rt: rate (Hz), integer preferred, typically, 1 .. 128, DEFAULT = 6;
%		Om: scale (cyc/oct), any real number, typically, .25 .. 4, DEFAULT = 2; 
%		Ph: (optional) symmetry (Pi) at f0, -1 < Ph < 1, DEFAULT = 0.
%	cond = (optional) [T0, f0, SF, BW, RO, df]; 
%		T0: duartion (sec), DEFAULT = 1.
%		f0: center freq. (Hz), DEFAULT = 1000.
%		SF: sample freq. (Hz), must be power of 2, DEFAULT = 16384
%		BW: excitation band width (oct), DEFAULT = 5.8.
%		RO: roll-off (dB/oct), 0 means log-spacing, DEFAULT = 0;
%		df: freq. spacing, in oct (RO=0) or in Hz (RO>0), DEFAULT = 1/16.
%	ph_c: component phase

% Acknowledge: This program is available due to Jian Lin's creative idea and
%	his C program [rip.c]. Thank Jonathan Simon and Didier Dipereux for their
%	Matlab program [ripfft.m].

% 08-Jun-98, v1.00

% default parameters and conditions
% para = [Am, Rt, Om, Ph]; cond = (optional) [T0, f0, SF, BW, RO, df];
para0 = [.9, 6, 2, 0]; cond0 = [1, 1000, 16384, 5.8, 0, 1/16];

% arguments
if nargin < 1, para = para0; end;
if nargin < 2, cond = cond0; end;
if nargin < 3, ph_c = []; end;
if length(para) < 4, para(4) = para0(4); end;
for k = 2:6, if length(cond) < k, cond(k) = cond0(k); end; end;

% parameter
Am = para(1); Rt = para(2); Om = para(3); Ph = para(4);

% excitation condition
f0 = cond(2);	% center freq
SF = cond(3);	% sample freq, 16384, must be an even number
BW = cond(4);	% bandwidth, # of octaves
RO = cond(5);	% roll-off at 3 dB/Oct, 0 means log-spacing
df = cond(6);	% freq. spacing, in oct (RO = 0) or in Hz (RO ~= 0)

% duration
T0 = cond(1); 	% actual duration in seconds
T1 = ceil(T0);	% duration for generating purpose
Ri = Rt*T1;	% modulation lag, # of df

% freq axis
if RO,	%compute tones freqs
	R1 = round(2.^([-1 1]*BW/2)*f0/df); fr = df*(R1(1):R1(2))';
else,	%compute log-spaced tones freqs
	R1 = round(BW/2/df); fr = f0*2.^((-R1:R1)*df)';
end;
M = length(fr);	% # of component
S = zeros(T1*SF/2, 1);	% memory allocation
% fprintf(2, 'Freq.: %d .. %d Hz\n', round(min(fr)), round(max(fr)));

fdx = round(fr*T1)+1;	% freq. index
x = log2(fr/f0);	% tono. axis (oct)

% roll-off and phase relation
r = 10.^(-x*RO/20);	% roll-off, -x*RO = 20log10(r)
if (~isempty(ph_c))
	th = ph_c;
else
	th = 2*pi*rand(M,1);	% component phase, theta
	disp('randomize component phase');
end
ph_c = th;
ph = (2*Om*x+Ph)*pi;	% ripple phase, phi, modulation phase

% modulation
S = [S; 0];
S(max(1, fdx-Ri)) = S(max(1, fdx-Ri)) + r.*exp(j*(th-ph));	% lower side
S(min(T1*SF/2+1, fdx+Ri)) = S(min(T1*SF/2+1, fdx+Ri)) + r.*exp(j*(th+ph));	% upper side
S = S * Am/2;
S(1) = 0; S = S(1:T1*SF/2);

% original stationary spectrum
S(fdx) = S(fdx) + r.*exp(j*th);  % moved here to save computation

% time waveform
s = ifft([S; 0; conj(flipud(S(2:T1*SF/2)))]);	% make it double side
s = real(s(1:round(T0*SF)));% only real part is good. *2 was ignored

