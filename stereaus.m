function s = stereaus(x, paras, NS)
% STEREAUS stereausis using IIR auditory filterbank (for band 180 -7246 Hz)
%	s = stereaus(x, [frmlen, tc, fac, shft]);
%	x	: the acoustic input, L_x-by-2 matrix.
%
%	COCHBA  = (global) [cochead; cochfil]; 
%		cochead : 1-by-M filter length (<= L) vector.
%		f  = real(cochead); filter order
%		CF = imag(cochead); characteristic frequency
%		cochfil : (Pmax+2)-by-M (L-by-M) [M]-channel filterbank matrix.
%		B = real(cochfil); MA (Moving Average) coefficients.
%		A = imag(cochfil); AR (AutoRegressive) coefficients.
%	M	: highest (frequency) channel 
%
%	PARAS	= [frmlen, tc, fac, shft];
%	frmlen	: frame length, typically, 8, 16 or 2^[natural #] ms.
%	tc	: time const., typically, 4, 16, or 64 ms, etc.
%		  if tc == 0, the leaky integration turns to short-term avg.
%	fac	: nonlinear factor (critical level ratio), typically, .01.
%		  The less the value, the more the compression.
%		  fac = 0,  y = (x > 0),   full compression, booleaner.
%		  fac = -1, y = max(x, 0), half-wave rectifier
%		  fac = -2, y = x,		 linear function
%	shft: shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
%		  etc. SF = 16K * 2^[shft].%	
%
%	STEREAUS computes the auditory spectrogram for an acoustic waveform.
%	This function takes the advantage of IIR filter's fast performance
%	which not only reduces the computaion but also saves remarkable
%	memory space.
%	See also: WAV2AUD, UNITSEQ

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
MASK = 0;

if MASK, mask = [0 0 0 0 1 0 0; 0 0 0 -3 0 1 0; 0 0 5 0 -3 0 1; ...
0 -3 0 5 0 -3 0; 1 0 -3 0 5 0 0; 0 1 0 -3 0 0 0; 0 0 1 0 0 0 0]; end;

% get filter bank,
%	L: filter coefficient length;
%	M: no. of channels
global COCHBA;
global TICKLABEL;
if nargin < 1,
	x = [cos(2*pi*(1:4000)'/8000*1000)+cos(2*pi*(1:4000)'/8000*800) ...
		cos(2*pi*(1:4000)'/8000*1010)+cos(2*pi*(1:4000)'/8000*810)];
end;
if nargin < 2, paras = [8, 0, -1, -1]; end;
if nargin < 3, NS = 1; end;

[L, M]	= size(COCHBA);	% p_max = L - 2;
L_x	= length(x);	% length of inpxt

% PARAS octave shift, nonlinear factor, frame length, leaky integration
shft	= paras(4);			% octave shift
%fac	= paras(3);			% nonlinear factor
fac	= -1;				% HWR
L_frm	= round(paras(1) * 2^(4+shft));	% frame length (points)

%paras(2) = 0;
if paras(2),
	alph	= exp(-1/(paras(2)*2^(4+shft)));	% decaying factor
else,
	alph	= 0;					% short-term avg.
end;

% get data, zero padding 
N = ceil(L_x / L_frm);		% # of frames
x(N*L_frm, 2) = 0;		% zero-padding
if nargout,
	s = zeros((M-1)*N, M-1);
else,
	s = [];
end;

% response and carry
y1 = zeros(L_frm, M-1);		% cochlear response
y2 = y1;
cy1 = zeros(L-2, M-1);		% carry (initial and final conditions)
cy2 = cy1;
cc = rand(M-1, M-1);	% cross channel correlation matrix (128 ch)

%CF = 440 * 2 .^ ((-31:97)/24);

% graphic stuff
% log. frequency axis
y_str = [];
for fdx = 5:-1:1,
		ystr = num2str(round(1000*2^(fdx-3+paras(4))));
		L1 = length(ystr);
		L2 = size(y_str, 2);
		if L1 < L2, ystr = [32*ones(1, L2-L1) ystr]; end;
		y_str = [ystr; y_str];
end;

if isa1map,
	h_cc = image(2:M, 2:M, real_col(cc));
else,
	noc = size(colormap, 1);
	h_cc = image(2:M, 2:M, cc/max(cc(:)));
end;
set(gca, 'ytick', 11:24:128, ['y' TICKLABEL], y_str, ...
		'xtick', 11:24:128, ['x' TICKLABEL], y_str);
ylabel('Channel #1: Frequency (Hz)');
xlabel('Channel #2: Frequency (Hz)');
grid on;
hold on;
plot([2 M], [2 M], 'k-.');
h_diag = plot(2:M, diag(cc)/max(cc(:)), 'b-');
h_disp = plot((-24:24)+35, zeros(49, 1), 'r-');
set([h_diag h_disp], 'linew', 3);
hold off;
axis xy;

for n = 1:N,
	
	x1 = x((1:L_frm)+(n-1)*L_frm, 1);
	x2 = x((1:L_frm)+(n-1)*L_frm, 2);

	for m = 1:M-1,
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ANALYSIS: cochlear filterbank
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% (IIR) filter bank convolution ---> y1
		p  = real(COCHBA(1, m+1));	% order of ARMA filter
		B  = real(COCHBA((0:p)+2, m+1));% moving average coefficients
		A  = imag(COCHBA((0:p)+2, m+1));% autoregressive coefficients
		[y1(:, m), cy1(1:p, m)] = filter(B, A, x1, cy1(1:p, m));
		[y2(:, m), cy2(1:p, m)] = filter(B, A, x2, cy2(1:p, m));

	end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% TRANSDUCTION: hair cells
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fluid cillia coupling (preemphasis) (ignored)

	% ionic channels (sigmoid function)
	y1 = sigmoid(y1, fac);
	y2 = sigmoid(y2, fac);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CORRELATION: cross-fiber
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% temporal integration window ---> y5
	if alph,	% leaky integration
		for k = 1:L_frm,
			cc = cc + alph * (y1(k, :)' * y2(k, :));
		end;
		cc = max(cc, 0);
	else,		% short-term average
		cc = 0;
		for k = 1:L_frm,
			cc = cc + (y1(k, :)' * y2(k, :));
		end;
		cc = cc / L_frm;
	end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DISPAIRITY PLOT 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if MASK,
		cc = max(filter2(mask, cc), 0);
	end;

	if ~rem(n, NS),
	if isa1map,
		set(h_cc, 'Cdata', real_col(cc));
	else,
	   	set(h_cc, 'Cdata', cc/max(cc(:))*noc);
		%set(h_cc, 'Cdata', cc, [0 max(cc(:))]);
	end;
	maxcc = max(cc(:));
	for k = -24:24,
		disp_cc(k+25) = sum(diag(cc, k));
	end;
	set(h_diag, 'Ydata', diag(cc)/maxcc*36+11);
	set(h_disp, 'Ydata', disp_cc/maxcc*2+84);
	drawnow;
	end;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% OUTPUT
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if nargout,
		s((1:M-1)+(M-1)*(n-1), :) = cc;
	end;

end;

