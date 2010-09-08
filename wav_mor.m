function [xmor, ymor] = wav_mor(x1, x2, para1, rv, sv)
% WAV_MOR cortical morphing (waveform-to-waveform version)
%	xmor = wav_mor(x1, x2);
%   [xmor, ymor] = wav_mor(x1, x2, para1, rv, sv);
%	x1, x2	: input, 2-column-matrix
%	xmor	: output, STEPS-column-matrix
%   para1 	= [paras STEPS ITER NORM FULLT FULLX]
%   paras	: [frame (ms), Tc (ms), Nonlinearity, shift (octave)], see WAV2AUD
%	STEPS   : # of morphing steps
%	ITER	: # of iteration for each morphing.
%   NORM	: normalization style, see COR2AUD, CORNORM
%	FULLT(X): overlapped output option, see COR2AUD
%	rv, sv  : rate, scale vector, see AUD2COR.
%
%   WAV_MOR morphs two sounds according some rules defined in CPLXMEAN.
%   See also: WAV2AUD, AUD2COR, COR2AUD, CPLXMEAN

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 15-Aug-97
% v1.01: 20-Aug-97, add non-truncation (FULL)
% v1.02: 12-Apr-98, remove non-causal option

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, bug: N = size(y, 2) --> size(y, 1)

% set default parameters
if nargin < 3, para1 = [8 4 -2 -1 9 30 1 0 0]; end; % 5 steps, 30 iteration
if nargin < 4, rv = 2 .^ (1:5); end;
if nargin < 5, sv = 2 .^ (-2:3); end;

% data length
% on/off set ramping (optional)
LH = 128;			% length of hanning window
HAN = .5*(1-cos(2*pi*(1:LH*2)'/(LH*2+1)));	% hanning window
L=length(x1); hdx=[1:LH (1:LH)+L-LH]; x1(hdx)=x1(hdx).*HAN; x1=unitseq(x1);
L=length(x2); hdx=[1:LH (1:LH)+L-LH]; x2(hdx)=x2(hdx).*HAN; x2=unitseq(x2);

% parameters adjustment
paras = para1(1:4);
if length(para1) < 5, STEPS = 5;	else, STEPS = para1(5); end;
if length(para1) < 6, ITER  = 30;	else, ITER  = para1(6); end;
if length(para1) < 7, NORM  = 1;	else, NORM  = para1(7); end;
if length(para1) < 8, FULLT = 0;	else, FULLT = para1(8); end;
if length(para1) < 9, FULLX = FULLT;else, FULLX = para1(9); end;

STF = 1000 / paras(1);		% sample temporal frequency
T = L/(16000*2^paras(4));	% duration
rv(find(rv<=1/T)) = [];		% adjust minimum rate
rv(find(rv>=STF/4)) = [];	% adjust maximum rate

% file name acq: tmp????1.cor and tmp????2.cor
fcor = 1;
while fcor > 0,
	fname(1, :) = ['tmp' setstr('a'+fix(rand(1, 4)*25)) '1.cor'];
	fcor = fopen(fname(1, :));
	if fcor > 0, fclose(fcor); end;
end;
fname(2, :) = [fname(1, 1:7) '2.cor'];	% for the second file
disp('IN CASE OF INCOMPLETE, YOU MAY LOAD ');
disp([cd '/' fname(1, 1:7)]); 
disp('FOR PARTIAL RESULTS.');

% auditory spectrum, rate-scale processing
y = wav2aud(x1, paras); maxN = size(y, 1);
aud2cor(y, [paras FULLT FULLX], rv, sv, fname(1, :)); 
y = wav2aud(x2, paras); maxN = max(size(y, 1), maxN);
aud2cor(y, [paras FULLT FULLX], rv, sv, fname(2, :));

% allocate memory for morphed data
L_frm = round(paras(1) * 2^(4+paras(4)));
xmor  = zeros(maxN*L_frm, STEPS);

%% morphing
x0 = rand(maxN*L_frm, 1);

for midx = 1:STEPS,
	MORPH = (midx-1)/(STEPS-1); 

	% morphing and spectrum reconstruction
	y = cor_mor(fname(1, :), fname(2, :), NORM, MORPH);
	y = aud_fix(y);	% fix the complex spectrum

	% keep the middle one
	if midx == fix((STEPS+1)/2), ymor = y; end;

	% waveform reconstruction
	N = size(y, 1); L = N*L_frm;
    x0=unitseq(interpft(x0,L))/2+(1-MORPH)*interpft(x1,L)+MORPH*interpft(x2,L);
    x0 = unitseq(x0);
	[x0, xmin, errv] = aud2wav(y, x0, [paras ITER 0 0]);
	
	% store the data 
	xmor(1:length(xmin), midx) = xmin;
	eval(['save ' fname(1, 1:7) ' xmor;']);	% in case of .. .. 
end;

%% killing tmp files
delete(fname(1, :));
delete(fname(2, :));
delete([fname(1, 1:7) '.mat']);
