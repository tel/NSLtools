function v5 = wav2aud(x, paras, filt, VERB)
% WAV2AUD fast auditory spectrogramm (for band 180 - 7246 Hz)
%	v5 = wav2aud(x, [frmlen, tc, fac, shft], filt, VERB);
%	x	: the acoustic input.
%	v5	: the auditory spectrogram, N-by-(M-1) 
%
%	COCHBA  = (global) [cochead; cochfil]; (IIR filter)
%       cochead : 1-by-M filter length (<= L) vector.
%               f  = real(cochead); filter order
%               CF = imag(cochead); characteristic frequency
%	cochfil : (Pmax+2)-by-M (L-by-M) [M]-channel filterbank matrix.
%		B = real(cochfil); MA (Moving Average) coefficients.
%		A = imag(cochfil); AR (AutoRegressive) coefficients.
%	M	: highest (frequency) channel 
%
%	COCHBA  = [cochfil]; (IIR filter)
%	cochfil : (L-by-M) [M]-channel filterbank impulse responses.
%
%	PARAS	= [frmlen, tc, fac, shft];
%	frmlen	: frame length, typically, 8, 16 or 2^[natural #] ms.
%	tc	: time const., typically, 4, 16, or 64 ms, etc.
%		  if tc == 0, the leaky integration turns to short-term avg.
%	fac	: nonlinear factor (critical level ratio), typically, .1 for
%		  a unit sequence, e.g., X -- N(0, 1);
%		  The less the value, the more the compression.
%		  fac = 0,  y = (x > 0),   full compression, booleaner.
%		  fac = -1, y = max(x, 0), half-wave rectifier
%		  fac = -2, y = x,         linear function
%	shft	: shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
%		  etc. SF = 16K * 2^[shft].%	
%
%	filt	: filter type, 'p'--> Powen's IIR filter (default)
%			       'p_o' --> Powen's old IIR filter (steeper group delay)	
%	
%	IIR filter : (24 channels/oct)
%	for the output of 	downsamp/shift	tc (64 ms)/ frame (16 ms)
%	==================================================================
%	180 - 7246		1	/0	1024	/ 256
%	90  - 3623		2	/-1	512	/ 128	*
%
%	Characteristic Frequency: CF = 440 * 2 .^ ((-31:97)/24);
%	Roughly, CF(60) = 1 (.5) kHz for 16 (8) kHz.
%
%	VERB	: verbose mode
%
%	WAV2AUD computes the auditory spectrogram for an acoustic waveform.
%	This function takes the advantage of IIR filter's fast performance
%	which not only reduces the computaion but also saves remarkable
%	memory space.
%	See also: AUD2WAV, UNITSEQ

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 04-Sep-98, add Kuansan's filter (as FIR filter)

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v2.00: 24-Jul-01, add hair cell membrane (lowpass) function

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v2.10: 04-Apr-04, remove FIR filtering option (see wav2aud_fir.m) 

% get filter bank,
%	L: filter coefficient length;
%	M: no. of channels

if nargin < 4, VERB = 0; end;
if nargin < 3, filt = 'p'; end;

if (filt=='k')
   error('Please use wav2aud_fir function for FIR filtering!');
end

if (filt == 'p_o') load aud24_old;
else global COCHBA; end;

[L, M] = size(COCHBA);	% p_max = L - 2;
L_x = length(x);	% length of input

% octave shift, nonlinear factor, frame length, leaky integration
shft	= paras(4);			% octave shift
fac	= paras(3);			% nonlinear factor
L_frm	= round(paras(1) * 2^(4+shft));	% frame length (points)

if paras(2),
	alph	= exp(-1/(paras(2)*2^(4+shft)));	% decaying factor
else,
	alph	= 0;					% short-term avg.
end;

% hair cell time constant in ms
haircell_tc = 0.5;
beta = exp(-1/(haircell_tc*2^(4+shft)));

% get data, allocate memory for ouput 
N = ceil(L_x / L_frm);		% # of frames
x(N * L_frm) = 0;		% zero-padding
x = x(:);
v5 = zeros(N, M-1);
%CF = 440 * 2 .^ ((-31:97)/24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% last channel (highest frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p	= real(COCHBA(1, M));
B	= real(COCHBA((0:p)+2, M));
A	= imag(COCHBA((0:p)+2, M)); 
y1	= filter(B, A, x); 
y2	= sigmoid(y1, fac);
% hair cell membrane (low-pass <= 4 kHz); ignored for LINEAR ionic channels
if (fac ~= -2), y2 = filter(1, [1 -beta], y2); end;
y2_h = y2;
y3_h = 0;

t0 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All other channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ch = (M-1):-1:1,

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ANALYSIS: cochlear filterbank
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% (IIR) filter bank convolution ---> y1
	p  = real(COCHBA(1, ch));	% order of ARMA filter
	B  = real(COCHBA((0:p)+2, ch));	% moving average coefficients
	A  = imag(COCHBA((0:p)+2, ch));	% autoregressive coefficients
	y1 = filter(B, A, x); 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% TRANSDUCTION: hair cells
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fluid cillia coupling (preemphasis) (ignored)

	% ionic channels (sigmoid function)
	y2 = sigmoid(y1, fac);

	% hair cell membrane (low-pass <= 4 kHz) ---> y2 (ignored for linear)
	if (fac ~= -2), y2 = filter(1, [1 -beta], y2); end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% REDUCTION: lateral inhibitory network
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% masked by higher (frequency) spatial response
	y3   = y2 - y2_h;
	y2_h = y2;

	% spatial smoother ---> y3 (ignored)
	%y3s = y3 + y3_h;
	%y3_h = y3;

	% half-wave rectifier ---> y4
	y4 = max(y3, 0);

	% temporal integration window ---> y5
	if alph,	% leaky integration
		y5 = filter(1, [1 -alph], y4);
		v5(:, ch) = y5(L_frm*(1:N));
	else,		% short-term average
        if (L_frm == 1)
    	    v5(:, ch) = y4;
        else
        	v5(:, ch) = mean(reshape(y4, L_frm, N))';
        end
	end;

	if (VERB & filt == 'p') , if rem(ch, 24) == 0,
		if exist('octband'), octband = octband + 1; 
		else, octband = 1; end;
		fprintf(2, '%d octave(s) processed\r', octband);
	end; end;

end;

if VERB, fprintf(2, '\n'); end;
