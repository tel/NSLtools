function y = aud_warp(y, sv, t21, g, distp, SRF, FULLOUT)
% AUD_WARP direct auditory mapping by cortical "warping" and "distance"
%	y = aud_warp(y, sv, t21, g1, distp);
%	y = aud_maps(y, sv, t21, g1, distp, SRF, FULLOUT);	
%	y	: source (mapped) auditory spectrogram (N-by-M)
%	sv	: scale vector (K2-by-1)
%	t21, g1	: magnitude warping path and gain (M-by-K2)
%   distp	: phase distance (M-by-K2)
%	SRF	: (optional) sample ripple frequency
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
if nargin < 6, CLIP = 0; end;   % default: difference, no clipping
if nargin < 7, FULLOUT = 0; end;% nonoverlapped output

[N, M]	= size(y);
LAP	= 20;
t0	= clock;

FULL = length(dista)/M - 1;

for n = 1:N,
	z = aud2cors(y(n, :), sv, SRF, FULL);		% cortical represent.
	z = cor_maps(z, t21, distp, g);			% distance mapping
	z = cor2auds(z, sv, 1, SRF, M, FULLOUT);	% auditory spectrum
	%y(n, :) = z(:)' - mean(z) + mean(y(n, :));	% replace mean
	y(n, :) = z(:).';
	time_est(n, N, LAP, t0);
end;
