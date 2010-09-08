function z = cplxmean(z1, z2, par, UNWRAP)
% CPLXMEAN complex mean
%	z = cplxmean(z1, z2, par);
%	z = cplxmean(z1, z2, par, UNWRAP);
%	z1, z2	: input complex matrices
%	par	: 0 < partio < 1 of z2
%	z       : complex "mean" matrix, e.g., (1-par)*z1+par*z2
%	UNWRAP	: (optional) unwrap phase or not
%
%	CPLXMEAN computes the partial "mean" (interpolation) for two 
%	complex matrices. COR_MOR uses this function to do the spatio-
%	temporal morphing. The "mean" here is not the arithematic mean.
%	More likely, it is the geometric mean or some kind of variation. 
%	Basically, the partio is constant all the way. However, if a negative
%	number is given, gradual partio will be applied.
%	The magnitude interpolation could be linear, logrithmic or sinusoidal.
%	The phase interpolation could be unwrap phase or shortest angle.
%	Since it is very hard to finalize the definition, this program
%	will be modified from time to time.
%	See also: COR_MOR

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 20-Jun-97, include phase unwrap option

% dimensions
[N, M] = size(z1);
if ~(size(z2, 1) == N & size(z2, 2) == M),
	error('Size mismatched');
end;

% options
if nargin < 3, par = .5; end;
if nargin < 4, UNWRAP = 0; end;

% magnitude "mean"
a1 = abs(z1); [dum, idx1] = max(a1);		% reference indices 1
a2 = abs(z2); [dum, idx2] = max(a2);		% reference indices 2

if par >= 0,	% constant morphing
	%z = cos(par*pi/2)*a1 + sin(par*pi/2)*a2;	% sinusoidal interp.
	%z = (1-par)*a1 + par*a2;			% linear interp.
	UNWRAP = 1;
	z = (a1.^(1-par)).*(a2.^par);	% geometric interp.
else,		% gradual morphing
	if N == 1,
		par1 = .5;
		%z = cos(par1*pi/2)*a1 + ...
                %                sin(par1*pi/2)*a2;
		%z = (a1 + a2) / 2;
		z = sqrt(a1 .* a2);
	else,
		R2 = -par;
		for n = 1:N,
			R1 = (n - 1) / (N - 1);
			R1 = (R1 - .5) * R2 * 2;
			R1 = sign(R1) * (1-exp(-abs(R1)))/(1-exp(-R2));
			par1(n) = .5 * (1 + R1);
			%z(n, :) = cos(par1(n)*pi/2)*a1(n, :) + ...
			% 	sin(par1(n)*pi/2)*a2(n, :);
			%z(n, :) = (1-par1(n))*a1(n, :) + par1(n)*a2(n, :);
			z(n, :) = exp((1-par1(n))*log(a1(n, :)) + ...
				 par1(n)*log(a2(n, :)));
		end;
	end;
end;

% phase "mean"
z1 = angle(z1);
z2 = angle(z2);

if UNWRAP,

	% global unwrap
	pi2 = 2 * pi;
	z1 = unwrap(z1);
	z2 = unwrap(z2);
	for m = 1:M,
		R1 = round(z1(idx1(m), m)/pi2);
		z1(:, m) = z1(:, m) - R1 * pi2;
		R1 = round(z2(idx2(m), m)/pi2);
		z2(:, m) = z2(:, m) - R1 * pi2;
	end;

	if par >= 0,
		z1 = (1-par) * z1 + par * z2;
	else,
		for n = 1:N,
			z1(n, :) = (1-par1(n))*z1(n, :) + ...
                                	par1(n)*z2(n, :);
		end;
	end;

else,

	% local unwrap
	z2 = exp(i*z2) ./ exp(i*z1);
	z2 = angle(z2);			% smallest phase difference

        if par >= 0,
                z1 = z1 + par * z2;
        else,
                for n = 1:N,
                        z1(n, :) = z1(n, :) + par1(n)*z2(n, :);
                end;
        end;

end;

% complex
z = z .* exp(i*z1);
