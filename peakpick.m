function y = peakpick(x, DIM);
% PEAKPICK peak picking.
%	y = peakpick(x);
%	y = peakpick(x, DIM); 
%
%	PEAKPICK picks the peaks of a sequence or a matrix.
%	If [x] is a matrix, y(i, :) indicates the peaks in i-th column for 1-D
%	or [y] indicates the peaks for 2-D (DIM == 2).

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 25-Jun-97

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 04-Sep-98, add window length parameter for 2D peak picking

if size(x, 1)==1, x = x'; end, 
[M, N] = size(x);
y = zeros(M, N);
if nargin < 2, DIM = 1; end;
wl = 2; % window length (1 sided) for determining the peak

if DIM == 1,
	x = sign(diff(x));
	x = sign(diff(x) + .5);
	y(2:M-1, :) = sign(1 - x);
else,
	x1 = zeros(M+2*wl, N+2*wl);
	for m = -wl:wl, for n = -wl:wl, if m | n,
		idx = (1:M) + m + wl;
		idy = (1:N) + n + wl;
		x1(idx, idy) = max(x1(idx, idy), x);
	end; end; end;
	y = x > x1([1:M]+wl, [1:N]+wl);
end;
