function [p, fr] = audpitch(y, paras);
% AUDPITCH pitch detection by observing auditogram
%	[p, fr] = audpitch(y, paras);
%	p:	pitch
%	fr:	frequency
%	y:	N-by-M auditory spectrogram
%	Contraints: 1. Trace upto fifth harmonic
%				2. Third harmonic has to be in the scope

% v1.00, 09-Sep-97

% template
pdx = [1:2, 25:26, 38:41, 48:51, 56:59];
pw  = [0.5560 0.4440 0.5464 0.4536 0.1846 0.3282 0.2798 0.2074 0.1955 ...
	0.3268 0.2744 0.2033 0.2439 0.3234 0.2442 0.1885];

[N, M] = size(y);
p = zeros(N, M);

K1 = [ones(1,12)*5 4 ones(1,23)*3 2 ones(1,91)];
KK = [ones(1,107)*16 15 14 13 ones(1,5)*12 11 10 9 ones(1,7)*8 7 6 5];

for m = 1:M,
	for k = K1(m):KK(m),
		R1 = m-38+pdx(k);
		if R1 > 128, m, k, end;
		p(:, m) = p(:, m) + pw(k) * y(:, R1);
	end;
end;

fr = 440 * 2 .^ (((1:128)-68)/24+paras(4));

