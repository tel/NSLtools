function sw = warping(ss, tt, g, lambda)
% WARPING 1-D warping interpolation (128-point)
%	sw = warping(s1, t21, g2, lambda);
%	sw = warping(ss, tt, gg, lambda);
%	s1: real positive, M-by-2 matrix
%	t21: warped axes
%	tt = t21+i*t12;
%	gg = g2+i*g1;
%	lambda: L-vector, 0 < lambda < 1.
%	Note: This function is not efficient when # > 256

% dimension
if  nargin < 4, lambda = .5; end;
[M, N] = size(ss);
L = length(lambda);
sm = zeros(M, L);

% axis warping
t11 = (1:M)';
t21 = real(tt);

if N == 2,
	g = max(ss);
	s2 = ss(:, 2) / g(1);
	s1 = ss(:, 1) / g(2);

	t12 = imag(tt);
    % two-way interpolation
    for l = 1:L,
        r = lambda(l);
        s21 = interp1(t11, s1, (1-r)*t11+r*t21);
		s12 = interp1(t11, s2, (1-r)*t12+r*t11);
        sw(:, l) = s12 * (1-r) * g(1) + s21 * r * g(2);
    end;

else,
	s1 = ss(:, 1);
	g1 = max(s1);
	s1 = s1 / g1;

	% one-way interpolation
	for l = 1:L,
		r = lambda(l);
		s2 = interp1(t11, s1, (1-r)*t11+r*t21);
		%sw(:, l) = s2;
		sw(:, l) = s2 * ((1-r) + r*g) * g1;
	end;
end;
