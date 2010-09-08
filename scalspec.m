function s = scalspec(z, S)
% SCALSPEC scale spectrum
%	s = scalspec(z);
%	s = scalspec(z, S);
%	z	: static cortical representation difference (M-by-K)
%	s	: scale spectrum (1-by-K)
%	S	: (optional) power option
%	
%	SCALSPEC return the squared scale response S for complex cortical 
%	response Z.
%	See also: AUD2CORS, AUD_PLTS

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 12-Nov-97, add power option

if nargin < 2, S = 2; end;
s = abs(z).^S;
s = sum(s);
s = (s).^(1/S);
