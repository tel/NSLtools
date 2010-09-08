function [H, pha] = minphase(mag)
% MINPHASE minimum phase response
%	[H, pha] = minphase(mag);
%	mag: (linear) single-side magnitde spectrum
%	pha: minimum phase for [mag].
%	H: minimum phase transfer function

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 25-Sep-97 

L = length(mag);

% avoid zero amplitude
if min(mag) <=0,
	mag1 = mag + 1e-5 * max(mag);
else,
	mag1 = mag;
end;

% log magnitude
R1 = log(mag1(:));
R1 = [R1; R1(L); R1(L:-1:2)];	% pure even magnitude function

% hilbert transform
%R1 = hilbert(R1);
R1 = fft(R1);
R1(2:L) = R1(2:L) * 2;		% select single side band
R1(L+1:2*L) = R1(L+1:2*L) * 0; 
R1 = ifft(R1);

% minimum phase
pha = -imag(R1);
pha = [0; (pha(2:L)-pha(2*L:-1:L+2))/2]; % pure odd phase function

% minimum phase transfer function
H	= mag(:) .* exp(i*pha);
