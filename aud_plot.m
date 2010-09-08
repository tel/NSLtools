function h = aud_plot(y, paras, Nor, NM)
% AUD_PLOT auditory spectrum plot 
%	aud_plot(y, paras);
%	aud_plot(y, paras, NM);
%	y: auditory spectrum (N-by-M)
%	paras: parameters (see WAV2AUD)
%	Nor: normalization value for the plot
%	NM: (optional) original dimension
%	h: image handle
%
%	AUD_PLOT plots the auditory spectrum Y in the current axis and
%	returns image handel H. The abscissa is the time axis in ms 
%	(linear scale) and the ordinate is the frequency axis in kHz
%	(logrithmic scale). Only two elements of PARAS are required, i.e.,
%	[frame length] = paras(1) and [octave shift] = paras(4). The
%	main purpose of this program is to use the same colormap as the
%	one which is used by the cortical representations.
%	See also: WAV2AUD, COR2AUD

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 13-Jun-97
% v1.01: 28-Jul-97, colormap auto-detection
% v1.02: 30-Jul-97, make it Matlab4 compatible
% v1.03: 24-Aug-97, plot an overlapped spectrogram

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 08-Sep-98, plot aud. spectrogram out of Kuansan's filter 
% v1.10: 21-Jun-99, display more than 5 oct in freq. axis
 
global TICKLABEL;

if ~isreal(y), y = aud_fix(y); end;

% dimensions
[N, M] = size(y);
if nargin < 4,
	FULLOUT = 0;
        ndx = (1:N) * paras(1);
        mdx = 1:M;
else,
	FULLOUT = 1;
	N0 = NM(1);
	M0 = NM(2);
	dN = fix((N - N0)/2);
	dM = fix((M - M0)/2);
	ndx = (-dN:N0+dN) * paras(1);
	mdx = -dM:M0+dM;
end;

% log. frequency axis
y_oct = floor(M/24);
y_str = [];
for fdx = y_oct:-1:1,
	ystr = num2str(round(1000*2^(fdx-(y_oct-2)+paras(4))));
	L1 = length(ystr);
	L2 = size(y_str, 2);
	if L1 < L2, ystr = [32*ones(1, L2-L1) ystr]; end;
	y_str = [ystr; y_str]; 
end;

% color mapping (while using A1MAP - 256 colors)
if isa1map,
	y = real_col(y); 
	h = image(ndx, mdx, y');	% no normalization
else,
	if nargin < 3, Nor = max(y(:)); end;
	h = imagesc(ndx, mdx, y', [0 Nor]);
end;
axis xy;

if FULLOUT, hold on;
	plot([0 N0+1 N0+1 0 0]*paras(1), ...
		[0 0 M0+1 M0+1 0], 'k--');
hold off; end;


% ylabel, xlabel
if (M == 95)
	y_str(1,:) = [];
	set(gca, 'ytick', 20:20:95, ['y' TICKLABEL], y_str);
else
	set(gca, 'ytick', 12:24:y_oct*24, ['y' TICKLABEL], y_str);
end
ylabel('Frequency (Hz)');
xlabel('Time (ms)');
