function [CF, B, A] = cochfil(n, shft);
% COCHFIL cochlear filter coefficient reader
%	[CF, B, A] = cochfil(n, shft);
%	HH = cochfil(CHAR);
%	n: channel indices, e.g., 60 or 11:100
%	shft: octave shift
%	CF: characteristic frequencies of the selected channels
%	B, A: IIR coefficients of the selected channel (only one)
%	HH: CHAR = 'o', overall response;
%		CHAR = 'n', normalized overall response.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00, 05-Sep-97
% v1.01, 01-Jan-98, added overall response

global COCHBA;

if nargin < 2, shft = 0; end;
if n(1) == 'o' | n(1) == 'n',	% overall response,
	[N, M] = size(COCHBA);
	CF = 440 * 2.^(((1:M)-31)/24)/16000;
	s = exp(2i*pi*CF);
	HH2 = 0;
	for m = 1:M,
		p = real(COCHBA(1, m));
		B = real(COCHBA(2+(0:p), m));
		A = imag(COCHBA(2+(0:p), m));
		h = polyval(B,s) ./ polyval(A,s);
		if n(1) == 'o',
			HH2 = HH2 + abs(h).^2;
		else,
			NORM = imag(COCHBA(1, m));
			HH2 = HH2 + abs(h).^2 / NORM;
		end;
	end;
	CF = HH2;

	B = [];
	A = [];
 

else,
	% filterbank
	if length(n) > 1,	% no coefficient output
		B = [];
		A = [];
	else,			% coefficient output
		if n > 0 & n < size(COCHBA, 2),
			p = real(COCHBA(1, n));
			B = real(COCHBA(2+(0:p), n));
			A = imag(COCHBA(2+(0:p), n));
		else,
			error('invalid channel index!');	
		end;
	end;

	% characteristic frequency
	CF = 440 * 2.^((n-31)/24+shft);
end;
