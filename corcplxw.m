function corcplxw(z, fcor)
% !CORCPLXW read complex numbers from a cortical file
%	corcplxw(z, fcor);
%	!!! NOT RECOMMENDED FOR EXTERNAL USE !!!
%	z: N-by-M complex matrix
%	fcor: the file handle
%
%	CORPLEXW writes N*M*2 numbers into a cortical
%	file.
%	See also: AUD2COR

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 11-Aug-97

% matrix size
[N, M] = size(z);
data_format = 'uint8';

% maximum absolute value
z = [real(z) imag(z)];
R1 = abs(z);
R1 = max(R1(:));
z  = lin2mu(z/R1);
fwrite(fcor, R1, 'float');

% write
count = fwrite(fcor, z, data_format);
if count < N*M*2, error('FileWritingError!'); end;
