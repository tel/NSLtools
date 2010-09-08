function [arg1, arg2] = corcplxr(fcor, N, M, MAGONLY)
% !CORCPLXR read complex numbers from a cortical file
%	z = corcplxr(fcor, N, M);
%	[head_size, data_size] = corcplxr;
%	!!! NOT RECOMMENDED FOR EXTERNAL USE !!!
%	z: N-by-M complex matrix
%	fcor: the file handle
%
%	CORPLEXR reads N*M*2 numbers in a cortical
%	file and then composes them into a N-by-M
%	complex matrix.
%	See also: AUD2COR

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 10-Jun-97
% v1.01: 11-Aug-97, mu-law decoding

if nargin < 4, MAGONLY = 0; end;
if nargout == 1, 
	data_format = 'uint8';
	R1 = fread(fcor, 1, 'float');
	if MAGONLY,
		z  = fseek(fcor, N(1)*M(1)*2, 'cof');
		arg1 = R1;
	else,
		z  = fread(fcor, [N(1) M(1)*2], data_format);
		z  = mu2lin(z) * R1;
		arg1 = z(:, 1:M(1)) + i*z(:, M(1)+1:M(1)*2);
	end;
	arg2 = [];
elseif nargout == 2,
	arg1 = 4;
	arg2 = 1;
end;
