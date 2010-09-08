function y1 = wav2y1(x, ch_sel, filt)
% WAV2Y1 fast auditory filter bank output
%	y1 = wav2y1(x, ch_sel, filt);
%	x	: the acoustic input.
%	y1	: the filter bank output, Lx-by-M 
%
%	COCHBA  = (global) [cochead; cochfil]; (IIR filter)
%		cochead : 1-by-M filter length (<= L) vector.
%			f  = real(cochead); filter order
%			CF = imag(cochead); characteristic frequency
%			M   : highest (frequency) channel
%		cochfil : (Pmax+2)-by-M (L-by-M) [M]-channel filterbank matrix.
%			B = real(cochfil); MA (Moving Average) coefficients.
%			A = imag(cochfil); AR (AutoRegressive) coefficients.
%
%	COCHBA  = [cochfil]; (FIR filter)
%	cochfil : (L-by-M) [M]-channel filterbank impulse responses.
%
%	shft	: shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
%		  etc. SF = 16K * 2^[shft].%	
%
%	filt	: filter type, 'p'--> Powen's IIR filter (default)
%		        	       'k'--> Kuansan's FIR filter
%	
%	IIR filter : (24 channels/oct)
%	for the output of 	downsamp/shift
%	===================================
%	180 - 7246			1	/0
%	90  - 3623			2	/-1	
%
%	Characteristic Frequency: CF = 440 * 2 .^ ((-31:97)/24);
%	Roughly, CF(60) = 1 (.5) kHz for 16 (8) kHz.
%
%	FIR filter : (20 channels/oct)
%	for the output of 	downsamp/shift	
%	===================================
%	258 - 6727	    	1	/0	
%	129 - 3363	    	2	/-1	
%
%	Characteristic Frequency: CF = 500 * 2 .^ ((-20:75)/20);
%	Roughly, CF(40) = 1 (.5) kHz for 16 (8) kHz.
%	
%	WAV2Y1 computes the auditory spectrogram for an acoustic waveform.
%	This function takes the advantage of IIR filter's fast performance
%	which not only reduces the computaion but also saves remarkable
%	memory space.
%	See also: AUD2WAV, UNITSEQ

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 02-Mar-98

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 08-Sep-98, add Kuansan's filter (as FIR filter)

% get filter bank,
%	L: filter coefficient length;
%	M: no. of channels

if nargin < 3, filt = 'p'; end;
if nargin < 2, ch_sel = -1; end;

if (filt == 'k') load aud20ks;
else global COCHBA; end;

[L, M] = size(COCHBA);	% p_max = L - 2;
L_x = length(x);	% length of inpxt


% channel selected
if ch_sel(1) < 0,	% increment
	ch_sel = 1:(-ch_sel):M;
else,
	ch_sel(find(ch_sel > M)) = [];
end;

% figure the output size
K = length(ch_sel);
N = L_x * K;
MB = 2^17;
if N > MB,
	disp(['Estimated Output Size: ' num2str(N/MB) 'MB']);
end;

x = x(:);
t0 = clock;
y1 = zeros(L_x, K);
LAP = max(1, min(K, ceil(MB/L_x)));
for k = 1:K,

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ANALYSIS: cochlear filterbank
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% (IIR) filter bank convolution ---> y1
	ch = ch_sel(k);
	if (filt == 'k')
	B	= COCHBA(:, ch);	% FIR filter bank
	A	= 1;
	else	
	p  = real(COCHBA(1, ch));	% order of ARMA filter
	B  = real(COCHBA((0:p)+2, ch));	% moving average coefficients
	A  = imag(COCHBA((0:p)+2, ch));	% autoregressive coefficients
	end
	y1(:, k) = filter(B, A, x); 
	time_est(k, K, LAP, t0);

end;
