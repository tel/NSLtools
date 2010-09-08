function y = aud_maps(y, sv, da, dp, SRF, dg, GRAD, FULLOUT)
% AUD_MAPS direct auditory mapping by cortical "distance"
%	y = aud_maps(y, sv, da, dp);
%	y = aud_maps(y, sv, da, dp, SRF, g, FULLOUT);	
%	y	: source (mapped) auditory spectrogram (N-by-M)
%	sv	: scale vector (K2-by-1)
%	da	: magnitude distance (M-by-K2)
%	dp	: phase dance (M-by-K2)
%	SRF	: (optional) sample ripple frequency
%	dg	: (optional)
%		length(g) > 1: gain vector
%		dg = 0: magnitude difference (default)
%		dg > 0: magnitude ratio, specifies lower clip level in dB
%	GRAD: gradual mapping
%	FULLOUT	: (optional) overlaped output
%
%	AUD_MAPS maps an auditory spectrogram into another according to 
%	the given cortical "distance" directly.
%	See also: WAV2AUD, COR_DIST

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 20-Aug-97, add non-truncation (FULL)
% v1.02: 05-Jan-98, bug: z=cor2auds(z,sv,.9,SRF*,M*,FULLOUT); and z(:)'

if nargin < 5, SRF = 24; end;	% 24 ch / oct
if nargin < 6, dg = 0; end;		% default: difference, no clipping
if nargin < 7, GRAD = 0; end;	% direct mapping
if nargin < 8, FULLOUT = 0; end;% nonoverlapped output

[N, M]	= size(y);
LAP	= 20;
t0	= clock;

%dM = (length(da)-M) /2;
FULL = length(da)/M - 1;

for n = 1:N,
	if GRAD,
		if n < N/3, lambda = 0;
		elseif n < N/3*2, lambda = (n-N/3)/(N/3); 
		else, lambda = 1; end;
	else,
		lambda = [];
	end;
	z = aud2cors(y(n, :), sv, SRF, FULL);		% cortical represent.
	z = cor_maps(z, da, dp, [dg lambda]);		% distance mapping
	z = cor2auds(z, sv, 1, SRF, M, FULLOUT);	% auditory spectrum
	%y(n, :) = z(:)' - mean(z) + mean(y(n, :));	% replace mean
	y(n, :) = z(:).';
	time_est(n, N, LAP, t0);
end;
