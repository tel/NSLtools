function zh = cor_maps(z1, da, dp, ARG4)
% COR_MAPS static mapping by cortical "distance"
%	zh = cor_maps(z1, da, dp, ARG4);
%	z1:	source complex matrix
%	da:	magnitude dance
%	dp:	phase dance
%	zh:	mapped complex matrix
%	ARG4: [g lambda] (optional)
%		length(g) > 1: gain vector
%		g = 0: magnitude difference (default)
%		g > 0: magnitude ratio, specifies lower clip level in dB
%		lambda: mapping coefficient, [0, 1].
%	COR_MAPS maps one frame of cortical representation to another.
%	See also: COR_MAP, COR_DIST

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
			 
if nargin < 4, CLIP = 0; end;		% default: difference, no clipping
% magnitude
a1 = abs(z1);

if length(ARG4) > 2,
	[M, K] = size(z1);
	if length(ARG4) > K,
		lambda = ARG4(K+1);
	else,
		lambda = 1;
	end;
	for k = 1:K,
		ah(:, k) = warping(a1(:, k), da(:, k), ARG4(k), lambda);
	end;
else,
	if length(ARG4) > 1,
		lambda = ARG4(2);
	else,
		lambda = 1;
	end;
	if ARG4(1),
		ah  = a1 .* da.^lambda;		% multiply magnitude ratio
	else,
		ah  = a1 + da * lambda;		% add magnitude difference
	end;
end;

% phase
zh = angle(z1) + dp * lambda;		% add phase difference

% new complex matrix
zh = ah .* exp(i * zh);

