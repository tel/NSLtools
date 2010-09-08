function y = schemat8(x, octshft, FCORNAME)
% SCHEMAT8  plot the schematic
%	y = schemat8(x);
%	y = schemat8(x, octshft, FCORNAME);
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

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Aug-97
% v1.01: 23-Sep-97, added octshft;
FUZZY = 0;

if nargin < 2, octshft = -1; end;
if nargin < 3, FCORNAME = 'nsltools.cor'; end;

% graphics setup 
loadload;
font_title = 10;

figsize([11 8.5]*.75);
clf;

% time, frequency
L = length(x);
paras(4) = octshft;
paras(1:2) = paras(1:2) / 2 ^(octshft+1);
SF = 16000 * 2^paras(4);
T = L / SF * 1000;

% waveform 
x = unitseq(x);
y = wav2aud(x, paras);

% rate-scale
hilo_r = ['Slow';'Fast'];
hilo_s = ['Coarse';'Fine  '];
%nepo = '-+';

sgn_sel = [-1 1];
%rv = [2 4 8 16 32] * 2^(octshft+1);
%sv = [.25 .5 1 2 4 8];
%r_sel = [2 4];	% 4, 16 Hz
%s_sel = [2 5];  % .5, 4 cyc/oct
rv = [4 16] * 2^(octshft+1);
sv = [.5 4];

WD = .35; 
HG = .12;

trng = 144;	 % ms
frng = 1.2;	   % oct
tres = 8;	   % ms
tlen = 64;	 % pts
ch_oct = 24;	% ch/oct
flen = 64;	 % ch

tdx = fix(trng/tres); tdx = (-tdx:tdx);
tdx0 = tdx * tres; tdx = tdx + tlen;
fdx = fix(frng*ch_oct); fdx = (-fdx:fdx);
fdx0 = fdx/ch_oct; fdx = fdx + flen;

% auditory spectrum to cortical representation
disp('GENERATING CORTICAL FILE');
aud2cor(y, [paras 0 0 0 1], rv, sv, FCORNAME);

for sgn_idx = 1:2,
for rdx = 1:2,
for sdx = 1:2,
	LF = .1 + (sdx-1) * (WD+.1);
	BT = .4 + (1.5-sgn_idx)*2*((rdx-1)*(HG+.03)+.1);

	% rate response
	R1 = gen_cort(rv(rdx), tlen, 1000/tres);
	r1 = fft(R1, tlen*2);
	r1 = r1([(1:tlen)+tlen 1:tlen]);

	% scale response
	R2 = gen_corf(sv(sdx), flen, ch_oct);
	r2 = fft(R2, flen*2);
	r2 = r2([(1:flen)+flen 1:flen]);

	% total response
	if sgn_idx == 1,
		r1 = conj(r1);
	end;
	r = r2 * r1.';
	r = r(fdx, tdx);

	%sub_c = axes('position', [LF, BT, WD, HG]);
	sub_c = axes('position', [LF+.2*WD, BT, WD*.8, HG]);
	z = cor_pick(FCORNAME, 3-sgn_idx, (rdx), (sdx));
	image(cplx_col(z, .2)');  axis xy
	set(sub_c, 'xtick', [], 'ytick', [], 'box', 'on');
	%sub_f = axes('position', [LF-WD*.025, BT+HG*.6, WD/4, HG/2]);
	sub_f = axes('position', [LF, BT, WD*.21, HG]);
    %image(real_col(real(r))); axis xy;
    image(cplx_col(r)); axis xy;
    [R1, R2] = size(r);
        text('position', [R2*0.1 R1*.9], 'str', ...
                [num2str(rv(rdx)) ' Hz'], 'fontwe', 'bo');
        text('position', [R2*0.1 R1*.75], 'str', ...
                [num2str(sv(sdx)) ' c/o'], 'fontwe', 'bo');

    set(sub_f, 'xtick', [], 'ytick', [], 'box', 'on');
	if sgn_idx == 1,
		sub_no = (2-rdx)*2+sdx;
	else,
		sub_no = (rdx-1)*2+sdx;
	end;
	if sub_no == 1,
		text('position', [-R2*.2 R1*1.2], 'str', setstr('A'+sgn_idx-1), ...
			'fontsi', 14, 'fontwe', 'bo');
	end;
		
	title(num2str(sub_no));
end;
end;
end;

if nargin < 2, delete(FCORNAME); end;

