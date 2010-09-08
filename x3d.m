function B = x3d(A, OVRLAP)
% X3D to convert a matrix to its 3-dim equivalent matrix
%	B = x3d(A);
%	B = x3d(A, OVRLAP);
%	A = [ a1 | a2 | .. | aN]
%	B = [ a1+1c | a2+2c | .. | aN+Nc]
%	   = A + [1 2 .. N)] * c
%	OVRLAP: overlap factor, the # of lines shared one space

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 02-Mar-96

[M, N] = size(A);
if nargin < 2,
	OVRLAP = 10;
end;

c = max(max(abs(A)))/(OVRLAP+.5);

B = A/c + ones(M,1)*(1:N); 
