function [s,profile] = multimvripfft(rippleList, cond);
% MULTIMVRIPFFT generates multiple moving ripples via FFT
%	s = multimvripfft();
%	s = multimvripfft(rippleList);
%	s = multimvripfft(rippleList, cond);
%   [s, profile] = multimvripfft(rippleList, cond);
%	rippleList = [Am1, w1, Om1, Ph1
%				 Am2, w2, Om2, Ph2
%				 ....
%      			 AmN, wN, OmN, PhN];
%		Am: modulation depth, 0 < Am < 1, DEFAULT = 1;
%		w: rate (Hz), integer preferred, typically, 1 .. 128, DEFAULT = 8;
%		Om: scale (cyc/oct), any real number, typically, .25 .. 4, DEFAULT = 1; 
%		Ph: (optional) sine-symmetry at f0 in radians, DEFAULT = 0.
%       N >= 1
%	cond = (optional) [T0, f0, BW, SF, CF, df, RO, AF, Mo]; 
%		T0: duartion (sec), DEFAULT = 1.
%		f0: lowest freq. (Hz), DEFAULT = 250.
%		BW: excitation band width (oct), DEFAULT = 5.
%		SF: sample freq. (Hz), should be power of 2, DEFAULT = 16384
%		CF: component-spacing flag, 1->log spacing, 0->harmonic, DEFAULT = 1
%		df: freq. spacing, in oct (CF=1) or in Hz (CF=0), DEFAULT = 1/16 or 1.
%		RO: roll-off (dB/oct), DEFAULT = 0 (CF=1) or 3 (CF=0)
%		AF: amplitude flag, 1->linear, 0->log (dB), DEFAULT = 1;
%		Mo: amp. total mod: 0<Mo<1 (AF=1); 0<Mo dB (AF=0); , DEF=0.9 or 10 dB
%		wM: Maximum temporal velocity to consider in Hz (DEFAULT = 120).
%	profile = spectro-temporal envelope *linear* profile; 
%
% By Powen Ru
% Acknowledge: This program is available due to Jian Lin's creative idea and
%	his C program [rip.c]. Thank Jonathan Simon and Didier Dipereux for their
%	Matlab program [ripfft.m].
%
% 08-Jun-98, v1.00
%
% 12/98 generalized by Jonathan Simon 
%    allows mulitple ripples
%    interface changed somewhat (f0 is lowest freq, Ph in radians, new defaults)
%    named multimvripfft
% 10/99 generalized by Tai-Shih & Jonathan Simon 
%    algorithm changed to compute envelope and take fft of that!
%    Almost as fast as Jian Lin's algorithm, but it generalizes to
%    passed envelopes and non-linear operations (e.g. exponentiation) on
%    envelopes.
%    Also allows non-integer ripple velocities, at the expense of requiring
%    that (w*T) be an integer for all w.
% 10/00 Fixed cases when pads < 0 and when w = Om = 0 . JZS

% default rippleList and conditions
% rippleList = [Am, w, Om, Ph]; cond = (optional) [T0, f0, BW, SF,   CF, df, RO, AF, Mo, wM];
rippleList0  = [1,  8, 1,  0];  cond0 =           [1, 125, 5, 16000, 1, 1/20, 0, 1, 0.9, 120];

% arguments
if nargin < 1, rippleList = rippleList0; end;
if nargin < 2, cond = cond0; end;
if size(rippleList,2) < 4, rippleList(:,4) = rippleList0(4); end;
if length(cond) >= 5
	if cond(5)==0, cond0(6) = 1; cond0(7) = 3; end
end
if length(cond) >= 8
	if cond(8)==0, cond0(9) = 10; end
end
for k = 2:10, if length(cond) < k, cond(k) = cond0(k); end; end;

% rippleList
Am = rippleList(:,1); w = rippleList(:,2); Om = rippleList(:,3); Ph = rippleList(:,4)-pi/2;

% excitation condition
T0 = cond(1); 	% actual duration in seconds
f0 = cond(2);	% lowest freq
BW = cond(3);	% bandwidth, # of octaves
SF = cond(4);	% sample freq, 16384, must be an even number
CF = cond(5);	% component-spacing flag, 1->log spacing, 0->harmonic
df = cond(6);	% freq. spacing, in oct (CF=1) or in Hz (CF=0)
RO = cond(7);	% roll-off in dB/Oct
AF = cond(8);	% amplitude flag, 1->linear, 0->log (dB)
Mo = cond(9);	% amp. total mod: 0<Mo<1 (Af=1); 0<Mo dB (Af=0)
wM = cond(10);	% amp. total mod: 0<Mo<1 (Af=1); 0<Mo dB (Af=0)

if abs(round(abs(w)*T0)-abs(w)*T0) > 1e-5
	error('Ripple Velocities not commensurate with Time')
end

max_rvel = max(max(abs(w)), 1/T0);
max_rfrq = max(max(abs(Om)), 1/BW);
t_step = 1/(16*max_rvel);
f_step = 1/(16*max_rfrq);
t_env_size = round((T0/t_step)/2)*2; % guarantee even number for fft components
f_env_size = round((BW/f_step)/2)*2; % guarantee even number for fft components
t_env = [0:t_env_size-1]  *t_step;
f_env = [0:f_env_size-1].'*f_step;

% Compute the maximum and the minimum of the envelope
profile	= zeros(f_env_size,t_env_size);
for row = 1:size(rippleList,1)
	rip_amp 	= Am(row);
	rip_vel		= w(row);
	rip_freq	= Om(row);
	rip_phase	= Ph(row)+pi/2;
	f_phase = 2*pi*rip_freq*f_env + rip_phase;
	t_phase = 2*pi*rip_vel*t_env;
	profile = profile + ...
		rip_amp*(sin(f_phase)*cos(t_phase)+cos(f_phase)*sin(t_phase));
	end
	
min_pro = min(min(profile));
max_pro = max(max(profile));
max_pro = max(max_pro, -min_pro);
if abs(max_pro)<1e-7; max_pro = sum(Am);end % if all w = Om = Ph = 0

profile = profile/max_pro;  % appropriately normalize

if AF==1
	t_step = 1/(4*max_rvel);% if linear, max_rvel is Nyquest freq for envelope, but 2 isn't enough.
else
	t_step = 1/(2*round(wM*T0)/T0);% otherwise, use wM as lowest freq for envelope, but 2
end
t_env_size = round((T0/t_step)/2)*2; % guarantee even number for fft components
t_env = [0:t_env_size-1]*t_step;
% freq axis
if CF==0	%compute harmonic tones freqs
	fr = df*(round(f0/df):round(2.^BW*f0/df)).';
else	%compute log-spaced tones freqs
	fr = f0*2.^((0:round(BW/df*2)/2-1)*df).';
end
f_env = log2(fr./f0);
f_env_size = length(fr);	% # of component

% Compute the maximum and the minimum of the envelope
profile	= zeros(f_env_size,t_env_size);
for row = 1:size(rippleList,1)
	rip_amp 	= Am(row);
	rip_vel		= w(row);
	rip_freq	= Om(row);
	rip_phase	= Ph(row)+pi/2;
	f_phase = 2*pi*rip_freq*f_env + rip_phase;
	t_phase = 2*pi*rip_vel*t_env;
	profile = profile + ...
		rip_amp*(sin(f_phase)*cos(t_phase)+cos(f_phase)*sin(t_phase));
end
profile = profile/max_pro;  % appropriately normalize
	
if AF==1
	profile = 1+profile*Mo; % shift so background = 1 & profile is envelope
else
	profile = 10.^(Mo/20*profile);
end

%%%%%%%%%%%%%%%%%%%%%%%% freq-domain AM %%%%%%%%%%%%%%%%%%%%%%%%%
L_sig = T0*SF;	% length of signal

% roll-off and phase relation
th = 2*pi*rand(f_env_size,1);	% component phase, theta
r = 10.^(-log2(fr/f0)*RO/20);		% roll-off, RO = 20log10(r)

S = zeros(1, L_sig);	% memory allocation
t_env_size_2 = t_env_size/2;
for m = 1:f_env_size
	f_ind = round(fr(m)*T0);
	S_tmpA = fftshift(fft(profile(m,:)))*exp(j*th(m))*r(m)/t_env_size*L_sig/2;
	padzerosonleft = f_ind - t_env_size_2 - 1;
	padzerosonright = L_sig/2 - f_ind - t_env_size_2;
	if ((padzerosonleft > 0) & (padzerosonright > 0) )
		S_tmpB = [zeros(1,padzerosonleft),S_tmpA,zeros(1,padzerosonright)];
	  elseif ((padzerosonleft <= 0) & (padzerosonright > 0) )
		S_tmpB = [S_tmpA(1 - padzerosonleft:end),zeros(1,padzerosonright)];
	  elseif ((padzerosonleft > 0) & (padzerosonright <= 0) )
		S_tmpB = [zeros(1,padzerosonleft),S_tmpA(1:end+padzerosonright)];
	end
	S_tmpC = [0, S_tmpB, 0, fliplr(conj(S_tmpB))];
	S = S + S_tmpC; % don't really have to do it all--know from padzeros which ones to do...
end

s = real(ifft(S));
s = s';
