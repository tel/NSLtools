function [rsw, rv, sv] = cor2rs(fname)
% COR2RS cortical rate-scale weighting coefficient
%	[rsw, rv, sv] = cor2rs(fname);
%	fname:		input .cor file
%	rsw = [sgnw(1)*rw*sw'; sgnw(2)*rw*sw'];
%
%	COR2RS returnss the ovelall rate, scale response. 
%	[rsw] will be a 2*K1-by-K2 matrix. 

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 06-Oct-97
% v1.01: 12-Apr-98, remove non-causal option

% read info 
fcor	= fopen(fname);
[para1, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor);

STF	= 1000 / para1(1);	% temporal sample frequency
SRF	= 24;			% spatial sample frequency
M(4)	= floor(M(1)/2*FULLX);  % dM
N(4)	= floor(N(1)/2*FULLT);  % dN

% rate-scale loop
for rdx = 1:K1,
	% rate filtering
	for sgn = [1 -1],
		for sdx = 1:K2,
			% load file (magnitude only)
			R1  = corcplxr(fcor, N(1)+2*N(4), M(1)+2*M(4), 1);
			% weighting
			%R1 = [real(R1) imag(R1)];
			%R1 = abs(R1); R1 = max(R1(:));
			rsw(rdx+(1-sgn)/2*K1, sdx) = R1;
		end;
	end;
end;

fclose(fcor);
