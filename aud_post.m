function y = aud_post(y)
% AUD_POST auditory post windowing
%	y = aud_post(y);
%	y: N-by-128, or 128-by-N
%	AUD_POST will remove the mean of auditory spectrum and then apply a
%	window on it (to remove the edge effect in the cortical representation.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 22-Dec-97 

% dimension
[N, M] = size(y);

if M == 128,
	TRANSPOSE = 0;
elseif N == 128,
	y = y.';
	TRANSPOSE = 1;
	[N, M] = size(y);
else,
	error('Not a auditory spectrogram');
end;

% windowing/weighting function
w = ones(1, M);
W = 6;
w([1:W (1:W)+128-W]) = .5*(1 - cos(2*pi*(1:2*W)/(2*W+1)));

% processing
for n = 1:N,
	y0 = mean(y(n, :));
	y(n, :) = y(n, :) - y0;	% mean removal
	y(n, :) = y(n, :) .* w + y0;			% windowing
end;

if TRANSPOSE, y = y.'; end;
