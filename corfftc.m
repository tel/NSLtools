function [Z_cum, HH] = corfftc(z, Z_cum, N, M, HR, HS, HH)
% !CORFFTC cortical fft and cumulation 
%	[Z_cum, HH] = corfftc(z, Z_cum, N, M, HR, HS, HH);
%	!!! NOT RECOMMENDED FOR EXTERNAL USE
%	z: N-by-M complex matrix
%	Z_cum: cumulated reverse response
%	HR, HS: rate, scale transfer function
%	HH: cumulated reverse filter transfer function
%	
%	CORFFTC is an internal routine to perform the 2-D reverse filtering
%	and cumulate the response which be used to reconstruct the auditory
%	spectrogram

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 11-Jun-97
% v1.01: 19-Aug-97, M(4) = dM, N(4) = dN
% v1.02: 30-Sep-97, add causal option
% v1.03: 12-Apr-98, remove non-causal option

% 2-D FFT
Z = zeros(N(3), M(2));

if M(4),
	z(:, [1:M(1)+M(4) (1:M(4))+M(3)-M(4)]) = ...
		z(:, [(1:M(1)+M(4))+M(4) 1:M(4)]);
else,
	z(1, M(3)) = 0;	% why? I forgot. Oh, it is zero padding
end;

for n = 1:N(1)+N(4)*2,
   	R1 = fft(z(n, :));
   	Z(n, :) = R1(1:M(2));
end;

for m = 1:M(2),
	Z(:, m) = fft(Z(:, m));
end;

% cumulation
R1 = HR*HS';
HH = HH + R1 .* conj(R1);
Z_cum = Z_cum + R1 .* Z;

