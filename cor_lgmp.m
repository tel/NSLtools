function [lag, map] = cor_lgmp(z1, z2, sv)
% COR_LGMP cortical "lag" and "mapping"
%	[lag, map] = cor_lgmp(z1, z2, sv);
%	z2 = cor_lgmp(z1, lag, map, partio);
%	z1: source complex matrix (M-by-K)
%	z2: target complex matrix (M-by-K)
%	lag: translation lag
%	map: direct mapping 
%	CLIP: (optional), 60 dB
%
%	COR_LGMP is a two-way function which, if the SECOND ARGUMENT is a matrix,
%	computes the "lag" and "mapping"; otherwise it will apply the "lag" and
%	"mapping"
%	between two static cortical representations.
%	sound to another by COR_MAP or COR_MAPS. If CLIP == 0, D_MAG
%	could be negative in some occassions.
%	See also: AUD2CORS, COR2AUDS, COR_MAP, COR_MAPS, AUD_MAPS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 05-Jun-98


if nargin < 3, CLIP = 0; end;		% default: difference, no clipping

% magnitude
a1 = abs(z1);
a2 = abs(z2);
g = max(a2) + i * max(a1);

if ~CLIP,
    % difference
    da  = a2 - a1;       % magnitude difference
elseif CLIP == -1,
	[M, K] = size(a1);
	da = zeros(M, K);
	t0 = clock;
	for k = 1:K,
		ss = [a1(:, k), a2(:, k)];
		tt = warppath(ss);
		da(:, k) = tt;
		time_est(k, K, 1, t0);
	end;
else,
    % ratio
    thr = 10 ^ (-CLIP/20);      % lower clip level
    a1  = max(a1, max(a1(:))*thr);  % clip magnitude 1
    a2  = max(a2, max(a2(:))*thr);  % clip magnitude 2
    da = a2 ./ a1;       % magnitude ratio
end;

% phase
dp = (unwrap(angle(z2)')-unwrap(angle(z1)'))';	% unwrapped phase difference
%z1 = exp(i*angle(z1));
%z2 = exp(i*angle(z2));
%dp = angle(z2 ./ z1);			% "smallest" phase difference

