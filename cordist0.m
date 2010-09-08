function D2 = cordist0(z1, z2, SUM, S)
% CORDIST0 Euclidean cortical "distance" (L2)
%	D2 = cordist2(z1, z2);
%	D2 = cordist2(z1, z2, SUM, S);
%
%	z1: complex cortical response matrix 1 (M-by-K)
%	z2: complex cortical response matrix 2 (M-by-K)
%	d:  Euclidean distance
%	SUM: (optional) sum up the distance. 
%
%	CORDIST2 computes the Eulicdean "distance" between two static cortical
%	representations.
%	See also: AUD2CORS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 16-Oct-97

if nargin < 3, SUM = 1; end;
if nargin < 4, S = 2; end;

[M, K] = size(z1);

% symmetry integration 
D2 = 0;
DR = 0;

ph = (-1:2) * pi / 2;
cp = cos(ph);
sp = sin(ph);
wp = [1 1 1 1];
wp = wp / sum(wp);

RR1 = 0;
for l = 1:length(ph),

	% real response
	R1 = real(z1)*cp(l) + imag(z1)*sp(l); 
	R2 = real(z2)*cp(l) + imag(z2)*sp(l);

	dR = R1 - R2;

	dR = dR - mean(dR(:));

	% distance cumulation	
	D2 = D2 + (abs(dR).^S) * wp(l);	% phase weighting
end;

D2 = D2 /K;

if SUM, D2 = sum(D2(:)); end;
D2 = D2 .^ (1/S);
