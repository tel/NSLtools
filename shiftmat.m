function p1 = shiftmat(p, idx, jdx, WRAP);
% p1 = shiftmat(p, idx, jdx);
% p1 = shiftmat(p, idx, jdx, WRAP);
% p: orignal matrix
% p1: shifted matrix
% idx, jdx: the amount of shift
% WRAP: (optional) wrap around the shifted entries 
% SHIFTMAT shifts a matrix by idx-row and jdx-column. The index can be any
% integer regardless of sign.

if nargin < 4, WRAP = 0; end;
[M, N] = size(p);

if WRAP,
	if idx,
		if idx > 0,
			p1 = p([M-(idx-1:-1:0) 1:M-idx], :);
		else,
			p1 = p([-idx+1:M 1:-idx], :);
		end;
	end;
	if jdx,
		if jdx > 0,
			p1 = p1(:, [N-(jdx-1:-1:0) 1:N-jdx]);
		else,
			p1 = p1(:, [-jdx+1:N 1:-jdx]);
		end;
	end;
else,
	p1  = 0 * p;
	p1(max(1+idx, 1):min(M, M+idx), max(1+jdx, 1):min(N, N+jdx)) = ...
		p(max(1, 1-idx):min(M-idx, M), max(1, 1-jdx):min(N-jdx, N));
end;
