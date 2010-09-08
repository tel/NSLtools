function [y, dil0] = aud_dil(y, dil, dil_style)
% AUD_DIL auditory spectrogram dilation
% 	y = aud_dil(y, dil);
%	[y, dil0] = aud_dil(y, dil, dil_style);
%	y: auditory spectrogram, N-by-M 
%	dil:  dilation factor, 0 < dil < Inf
%	dil_style: (optional) dilation center style
%		= 0, at the middle channel (default)
%		= 1, at the most-energy channel
%		= 2, at the highest peak of each temporal row
%	dil0: dilation center, N-by-1
%	
%	AUD_DIL dilates the auditory spectrogram by a factor of [dil]. 
%	This function operates row-wise and will keep the energy intact.
%	See also: WAV2AUD, AUD2WAV

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

[N, M] = size(y);

% dilation center
if nargin < 3, dil_style = 0; end;
if dil_style == 1,
	[dum, dil0] = max(sum(y.^2));
	dil0 = dil0 * ones(1, N);
elseif dil_style == 2,
	[dum, dil0] = max(y.');
	[dum, idx] = max(dum);
	% only one channel shift allowed
	for k = idx+1:N,
		dil0(k) = max(min(dil0(k), dil0(k-1)+1), dil0(k-1)-1);;
	end; 
	for k = idx-1:-1:1,
                dil0(k) = max(min(dil0(k), dil0(k+1)+1), dil0(k+1)-1);; 
        end; 
else,
	dil0 = round((M+1)/2) * ones(1, N); 
end;
		
% dilation 
t0 = clock;
for n = 1:N,
	yn = y(n, :);

	% energy
	g = sum(yn);

	% interpolation
	ndx1 = (1:M)-dil0(n);
	ndx2 = min(max(ndx1/dil, ndx1(1)), ndx1(M));
	yn = interp1(ndx1, yn, ndx2);

	% normalization
	y(n, :) = yn / max(yn) * g;
	time_est(n, N, 60, t0);

end;
