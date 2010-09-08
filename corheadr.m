function [paras, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor)
% !CORHEADR read cortical file header
%	[para1, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor);
%	!!! NOT RECOMMENDED FOR EXTERNAL USED !!!
%	fcor	: handle of the cortical file
%	para1	: parameter set, [L_frm, tc, fac, shft], see WAV2AUD
%	K1 (K2)	: the length of rate (scale) vector.
%	rv (sv)	: rate (scale) vector
%	N, M	: dimensions of the original auditory spectrum
%   FULLT (FULLX): fullness of temporal (spectral) margin. The value can
%		  be any real number within [0, 1]. If only one number was
%		  assigned, FULLT = FULLX will be set to the same value.
%
%	CORHEADR reads cortical header into variables. This function
%	serves as the internal routine in many cortical manipulation
%	routines.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 19-Aug-97, added FULL
% v1.02: 28-Aug-97, remove 2-D mean Yn, Ym
% v1.03: 29-Sep-97, causal option
% v1.04: 12-Apr-98, remove non-causal option

% parameter
paras = fread(fcor, 4, 'float');
K1	  = fread(fcor, 1, 'float');
K2	  = fread(fcor, 1, 'float');
rv	  = fread(fcor, K1,'float');
sv	  = fread(fcor, K2,'float');
N	  = fread(fcor, 1, 'float');
M	  = fread(fcor, 1, 'float');
FULLT = fread(fcor, 1, 'float');
FULLX = fread(fcor, 1, 'float');

% spatial, temporal zeros padding
N1 = 2^nextpow2(N);	 N2 = N1*2;
M1 = 2^nextpow2(M);	 M2 = M1*2;

% N, M
N = [N; N1; N2];	% column vector
M = [M; M1; M2];
