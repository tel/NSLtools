function [da, dp, g] = cor_dist(z1, z2, CLIP)
% COR_DIST cortical "distance"
%	[da, dp] = cor_dist(z1, z2);
%	[da, dp, g] = cor_dist(z1, z2, CLIP);
%	z1: source complex matrix (M-by-K)
%	z2: target complex matrix (M-by-K)
%	da: magnitude distance (diff/ratio)
%	dp: phase distance
%	CLIP: (optional)
%	 = 0: magnitude difference (default)
%	 > 0: magnitude ratio, CLIP specifies lower clip level (dB)
%	 = -1: magnitude nonlinear warping
%	COR_DIST computes the "distance" between two static cortical
%	representations which, presumably, can be used to "map" one
%	sound to another by COR_MAP or COR_MAPS. If CLIP == 0, D_MAG
%	could be negative in some occassions.
%	See also: AUD2CORS, COR2AUDS, COR_MAP, COR_MAPS, AUD_MAPS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 22-Jan-98, included magnitude nonlinear warping
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

