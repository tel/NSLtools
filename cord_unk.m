function D2 = cord_unk(z1, z2, W, S)
% CORD_UNK cortical "distance" (L2) related to unknown weighting function
%	D2 = cordist2(z1, z2);
%	d2 = cordist2(z1, z2, W, S);
%	z1: complex cortical response matrix 1 (M-by-K)
%	z2: complex cortical response matrix 2 (M-by-K)
%	W (matrix): weighting function. => d2 (scalar): distance 
%	  (scalar): # of phase channel [P=16] => D2 (P-by-K): distance matrix
%	CORD_UNK computes the Eulicdean "distance" between two static cortical
%	representations provided the weighting function is given. Otherwise,
%	
%	See also: AUD2CORS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 16-Oct-97

[M, K] = size(z1);
% butter, 2nd, 3.12
Bx = [0.2113    0.2113];
Ax = [1.0000   -0.5774];

if nargin < 3, W = 16; end;
if nargin < 4, S = 2; end;

if sum(size(W)) > 2, 
	P = size(W, 1);
else,
	P = W;
end;
D2 = zeros(P, K);

% phase channel
ph = (0:P-1)/P * 2 * pi;
cp = cos(ph); sp = sin(ph);

for p = 1:P,

	% real response
	R1 = real(z1)*cp(p) + imag(z1)*sp(p); 
	R2 = real(z2)*cp(p) + imag(z2)*sp(p);

	%R1 = max(R1, 0); R2 = max(R2, 0);	
	
	% difference
	dR = R1 - R2;

	%dR = filter(Bx, Ax, dR);

	% distance integration along tono
	D2(p, :) = sum(abs(dR).^S);
end;

if sum(size(W)) > 2,
	D2 = D2.*W; D2 = sum(D2); D2 = sum(D2); D2 = D2.^(1/S);
end;
