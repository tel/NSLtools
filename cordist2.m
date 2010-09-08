function D2 = cordist2(z1, z2, SUM, S)
% CORDIST2 Euclidean cortical "distance" (L2)
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

D2 = zeros(M, K);

Bx = [0.1207   -0.1685    0.2655   -0.1685    0.1207];
Ax = [1.0000   -2.4760    2.8920   -1.6546    0.4185];

ph = (1:18)/9*pi;
cp = cos(ph);
sp = sin(ph);
L = length(ph);

for l = 1:L,

	% real response
	R1 = real(z1)*cp(l) + imag(z1)*sp(l); 
	R2 = real(z2)*cp(l) + imag(z2)*sp(l);

	% HWR
	%R1 = max(R1, 0); R2 = max(R2, 0);

	% difference
	dR = R1 - R2;

	% balance
	%dR = dR - mean(dR(:));
	%dR = filter(Bx, Ax, dR);

	%for k = 1:K,
	%	dR(:, k) = dR(:, k) - mean(dR(:, k));
	%end;

	% distance cumulation	
	D2 = D2 + (abs(dR).^S);
end;
D2 = D2 / L;

if SUM,
	D2 = mean(D2');
	D2 = mean(D2).^(1/S);
end;
