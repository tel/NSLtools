function yh = cornorm(Z_cum, HH, N, M, NORM, FOUTT, FOUTX)
% !CORNORM cortical normalization and ifft
%	yh = cornorm(Z_cum, HH, N, M, NORM, FOUTT, FOUTX);
%	!!! NOT RECOMMENDED FOR EXTERNAL USE !!!
%	Z_cum	: cumulated reverse response
%	HH	: overall 2-D response
%	N, M	: dimensions of the original auditory spectrum
%	%Yn, Ym	: 2-D "dc"
%	NORM	: 0=flat, 1=full, .x=partial normalization
%   FOUTT(X): (optional) overlapped output, within [0, 1]
%	% CAU: Why there is not causal option is this stage ?
%	% ans: Even you use causal temporal filter, you will need
%	%	non-causal temporal filters (i.e., reverse version
%	%	of the analysis filters) to serve as the synthesis
%	%	filters to undergo reconstruction. Thus we assume
%	%	the output spread is symmetric.
%
%	CORNORM normalizes 2-D rate-scale response and then does
%	the 2-D inverse fft to conclude reconstruction.
%	See also: COR2AUD, COR_SEL, COR_MAP, COR_MOR

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 20-Aug-97, modified NORM, added FULLOUT
% v1.02: 28-Aug-97, removed 2-D mean Yn, Ym
% v1.03: 02-Oct-97, add causal option

if nargin < 6, FOUTT = 0; end;
if nargin < 7, FOUTX = 0; end;

% modify overall transfer function
sumH = sum(HH(:));
HH = NORM * HH + (1 - NORM) * max(HH(:));
HH = HH / sum(HH(:)) * sumH;

% normalization
ndx = 1:N(3);
mdx = 1:M(2);
Z_cum(ndx, mdx) = Z_cum(ndx, mdx) ./ HH(ndx, mdx);

N(4) = floor(N(1)/2*FOUTT);					% dN
ndx	 = 1:N(1)+2*N(4);						% 1:N(1), if N(4)=0.
ndx1 = [(1:N(4))+N(3)-N(4) 1:N(1)+N(4)];	% 1:N(1), if N(4)=0.
M(4) = floor(M(1)/2*FOUTX);					% dM
mdx1 = [(1:M(4))+M(3)-M(4) 1:M(1)+M(4)];	% 1:M(1), if M(4)=0.

y 	= zeros(N(1)+2*N(4), M(2));				% zeros(N(1), M(2));
yh	= zeros(N(1)+2*N(4), M(1)+2*M(4));		% zeros(N(1), M(1));

% 2-D IFFT
for m = 1:M(2),
	R1 = ifft(Z_cum(:, m));
	y(:, m) = R1(ndx1);
end;
for n = ndx,
	R1 = ifft(y(n, :), M(3));
	yh(n, :) = R1(mdx1);
end;

yh = yh * 2;

% fix auditory spectrum
%yh = aud_fix(yh);
