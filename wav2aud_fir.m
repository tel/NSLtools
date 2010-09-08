function v5 = wav2aud_fir(x, paras, filt, VERB)
% WAV2AUD_FIR auditory spectrogram by FIR filtering (for band 250 - 6727 Hz)
%	v5 = wav2aud_fir(x, [frmlen, tc, fac, shft], filt, VERB);
%	x	: the acoustic input.
%	v5	: the auditory spectrogram, N-by-(M-1) 
%
%	COCHBA  = [cochfil]; (FIR filter)
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
%	filt	: filter type, 'k'--> Kuansan's FIR filter
%	
%	FIR filter : (20 channels/oct)
%	for the output of 	downsamp/shift	tc (64 ms)/ frame (16 ms)
%	==================================================================
%	258 - 6727		1	/0	1024	/ 256
%	129 - 3363		2	/-1	512	/ 128	*
%
%	Characteristic Frequency: CF = 500 * 2 .^ ((-20:75)/20);
%	Roughly, CF(40) = 1 (.5) kHz for 16 (8) kHz.
%	
%	VERB	: verbose mode
%
%	WAV2AUD_FIR computes the auditory spectrogram for an acoustic
%	waveform by FIR filtering!
%	See also: WAV2AUD
%	Acknowledge: Thank Elena Grassi for suggesting the fftfilt
%		     function to reduce FIR filtering processing time

% get filter bank,
%	L: filter coefficient length;
%	M: no. of channels

if nargin < 4, VERB = 0; end;
if nargin < 3, filt = 'k'; end;

if (filt == 'k') load aud20ks;
else error('Please use wav2aud function for IIR filtering!');  end;

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
B	= COCHBA(:, M);
y1	= fftfilt(B, x); 
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
	B	= COCHBA(:, ch);	% FIR filter bank
	y1 = fftfilt(B, x); 
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

	if VERB , if rem(ch, 20) == 0,
		if exist('octband'), octband = octband + 1; 
		else, octband = 1; end;
		fprintf(2, '%d octave(s) processed\r', octband);
	end; end;

end;

if VERB, fprintf(2, '\n'); end;
