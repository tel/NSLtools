function schematf(x, octshft, FCORNAME, SMALL)
% SCHEMATF  plot the schematic (forwarding part)
%	xh = schematf(x);
%	xh = schematf(x, octshft, FCORNAME);
%	x: input at 8 kHz
%
%	This function present an overview of the NSLtools
%	for 1 second, 8000 Hz sample rate, 8 ms frame
%	rv = 2.^(1:5); sv = 2.^(-2:3);
%	=============================================
%	6200 seconds for Quadra-650 (Mac)
%	1180 seconds for Power PC (Mac)
%	430 seconds for SPARCS 5
%	360 seconds for SPARCS 20
%	200 seconds for ULTRA-1
%	120 seconds for ULTRA-30

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Aug-97
% v1.01: 23-Sep-97, added octshft;
if nargin < 1, x = auread('_come.au'); end;
if length(x) < 10, x = auread('_come.au'); end;
if nargin < 2, octshft = -1; end;
if nargin < 3, FCORNAME = 'nsltools.cor'; end;
if isempty(FCORNAME), FCORNAME = 'nsltools.cor'; end;
if nargin < 4, SMALL = 0; end;

% graphics setup 
loadload;
colormap(1-gray);
arrow_co = 'r';		% arrow color
arrow_lw = 1;		% arrow linewidth
if SMALL, fontsi = 6; fig_ratio = .375;
else, fontsi = 10; fig_ratio = .75; end;	% title font
figsize([11 8.5]*fig_ratio);	% set figure size
clf;

% time, frequency
L = length(x);
paras(4) = octshft;							% octave shift
paras(1:2) = paras(1:2) / 2 ^(octshft+1);	% frame size and step
paras(3) = -2;								% linear
SF = 16000 * 2^octshft;						% sample rate
T = L / SF * 1000;

% waveform 
sub_x = axes('position', [.07 .75 .4 .15]);
x = unitseq(x);
h = plot((1:L)/SF*1000, x); 
set(h, 'linew', .1);
axis([0 L/SF*1000 -9 9]);
%text('position', [-.08*T, 11], 'str', '(a)');
text('position', [-.08*T, 11], 'str', 'A');
set(gca, 'fontsi', round(fontsi*.6));
title('Time Waveform'); xlabel('Time (ms)');
set(get(gca, 'title'), 'fontsi', fontsi);

% arrow box
sub_a0 = axes('position', [.48 .75 .04 .15]);
h_a0 = plot([0 1], [.5 .5], arrow_co, ...
	 [.8 1 .8], [.475 .5 .525], arrow_co);
set(h_a0, 'linew', arrow_lw);
axis([0 1 0 1]); axis off;

% spectrum
sub_y = axes('position', [.53 .75 .4 .15]);
text('position', [.5 .5], 'str', 'Please Wait .. ..', ...
		'fontsi', fontsi*2, 'co', 'w', 'ho', 'ce');
axis off;
drawnow;

y = wav2aud(x, paras);
aud_plot(y, paras);
set(gca, 'fontsi', round(fontsi*.6), 'yaxisloc', 'right');
set(get(gca, 'ylabel'), 'fontsi', round(fontsi*.6));
set(get(gca, 'xlabel'), 'fontsi', round(fontsi*.6));
title('Auditory Spectrogram');
set(get(gca, 'title'), 'fontsi', fontsi);
%text('position', [-.08*T, 1.1*128], 'str', '(b)');
text('position', [-.08*T, 1.1*128], 'str', 'B');

% arrow box
xa = [.795 .535 .29 .03];
R1 = [-1; 0; 1] * .005; R1 = R1*ones(1, 4);
sub_a1 = axes('position', [.07 .65 .86 .08]);
h_a1 = plot([[1;1]*xa [xa(1); xa(4)]], ...
		[[1; 0] [.4; 0]*[1 1 1] [.4;.4]], arrow_co, ...
		R1+[1;1;1]*xa, .1*[1;0;1]*ones(1, 4), arrow_co);
set(h_a1, 'linew', arrow_lw);
axis([0 1 0 1]); axis off;
%text('position', [-.03 -.1], 'str', '(c)');
text('position', [-.03 -.1], 'str', 'C');
text('position', [.4 .62], 'str', ...
	'Multiresolution Cortical Filters and Outputs', ...
	'ho', 'ce', 'fontsi', round(fontsi*1.2), 'co', 'b');

% rate-scale
hilo_r = ['Slow';'Fast'];
hilo_s = ['Coarse';'Fine  '];
nepo = '-+';

sgn_sel = [-1 1];
rv = [4 16] * 2^(octshft+1);
sv = [.5 2];
r_sel = [1 2];	% 4, 16 Hz
s_sel = [1 2];  % .5, 4 cyc/oct

WD = .4; 
HG = .15;

trng = 256;	 % ms
frng = 1.2;	 % oct
tres = 8;	 % ms
tlen = 64;	 % pts
ch_oct = 24; % ch/oct
flen = 32;	 % ch

tdx = fix(trng/tres); tdx = (1:tdx);
tdx0 = tdx * tres;
fdx = fix(frng*ch_oct); fdx = (-fdx:fdx);
fdx0 = fdx/ch_oct; fdx = fdx + flen;

sub_n = axes('position', [.07 .475 .4 .05]);
text('position', [.5 .5], 'str', 'Please Wait .. ..', ...
		'fontsi', fontsi*2, 'co', 'w', 'ho', 'ce');
axis off;
drawnow;

% auditory spectrum to cortical representation
disp('GENERATING CORTICAL FILE');
aud2cor(y, [paras 0 0 1], rv, sv, FCORNAME);

cla;
text('position', [.5 -2.25], 'str', 'Upward Moving', ...
		'fontsi', round(fontsi*1.4), 'co', 'b', 'ho', 'ce');
sub_p = axes('position', [.53 .475 .4 .05]);
text('position', [.5 -2.25], 'str', 'Downward Moving', ...
		'fontsi', round(fontsi*1.4), 'co', 'b', 'ho', 'ce');
axis off;


for sgn_idx = 1:2,
for rdx = 1:2,
for sdx = 1:2,

	if rdx ~= sdx,
	LF = .07 + (sgn_idx-1)*.46;
	BT = .1 + (2-sdx) * .3;

	% rate response
	R1 = gen_cort(rv(r_sel(rdx)), tlen, 1000/tres);
	r1 = ifft(R1, tlen*2);
	r1 = r1(1:tlen);

	% scale response
	R2 = gen_corf(sv(s_sel(sdx)), flen, ch_oct);
	r2 = ifft(R2, flen*2);
	r2 = r2([(1:flen)+flen 1:flen]);

	% total response
	if sgn_idx == 1,
		r1 = conj(r1);
	end;
	r = r2 * r1.';
	%fdx, tdx
	r = r(fdx, tdx);

	sub_f = axes('position', [LF, BT+HG, WD/4, HG/2]);
	%image_pm(real(r), [], [16 16], 4); axis xy;
	image_q(r, [], [16 16]);	% .0001 -- .001  --.1
	[R1, R2] = size(r);
	text('position', [R2*1.3 R1*3/4], 'str', ...
				[hilo_r(rdx, :) ' Rate'], 'fontsi', fontsi);
	text('position', [R2*1.3 R1/3], 'str', ...
				[deblank(hilo_s(sdx, :)) ' Scale'], 'fontsi', fontsi);

	set(gca, 'xtick', [], 'ytick', [], 'box', 'on');

	sub_c = axes('position', [LF, BT, WD, HG]);
	z = cor_pick(FCORNAME, 3-sgn_idx, r_sel(rdx), s_sel(sdx));

	%image_c(z.', mean(y(:))/4, [20 48], 4);
	image_q(z.', mean(y(:))/4, [20 48]);
	axis xy
	set(gca, 'xtick', [], 'ytick', [], 'box', 'on');
	end;	
end;
end;
end;
