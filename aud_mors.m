function [y, z1, z2, R1] = aud_mors(y1, y2, N, sv, METHOD)
% AUD_MORS direct auditory morphing
%	y = aud_mors(y1, y2);
%	y = aud_mors(y1, y2, N, sv, METHOD);	
%	y1, y2	: source auditory spectra (N-by-M)
%	N: # of morphing steps
%	sv	: scale vector (K2-by-1)
%	METHOD: 'warping' or 'mappning'
%
%	AUD_MORS morphs two auditory spectra into one spectrogram 
%	See also: WAV2AUD, AUD2CORS, COR2AUDS

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 27-Jan-98

if nargin < 3, STEPS = 5; end;	
if nargin < 4, sv = 2.^(-2:3); end;
disp('Looking for warping path .. ..')
M = length(y1);
K = length(sv);
SRF = 24;
FULL = .5;
M0 = 128;

z1 = aud2cors(y1, sv, SRF, FULL);	a1 = abs(z1);
z2 = aud2cors(y2, sv, SRF, FULL);	a2 = abs(z2);
if METHOD(1) == 'w', tt = cor_dist(z1, z2, -1); end;
p1 = unwrap(angle(z1)')';
p2 = unwrap(angle(z2)')';

t0 = clock;
LAP = 10;

disp('Morphing or Mapping .. ..');

for n = 1:N,
	lambda = n / (N+1);

	if METHOD(1)=='w',
		for k = 1:K,
			R1(:, k) = warping([a1(:, k) a2(:, k)], tt(:, k), [], lambda);
		end;
	else,
		R1 = a1.^(1-lambda) .* a2.^lambda;
	end;
	R1 = R1 .* exp(i*(p1*(1-lambda)+p2*lambda));
	yh = cor2auds(R1, sv, 1, SRF, M0);	% auditory spectrum
	y(n, :) = yh(:).';
	time_est(n, N, LAP, t0);
end;
