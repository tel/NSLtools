function cog = cors_cog(z)
% CORS_COG Euclidean cortical "center of gravity" (L2)
%	xpyc0 = cors_cog(z);
%
%	z: complex cortical response matrix (M-by-K)
%	xpyc: centroid 
%
%	CORS_0 computes the "centroid" of a static cortical representations.
%	See also: AUD2CORS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 24-Feb-98

[M, K] = size(z);
a = abs(z);
a2 = a.^2;
r_sum_a2 = sum(a2');
c_sum_a2 = sum(a2);
sum_a2 = sum(c_sum_a2);

if sum_a2,
	x0 = (1:M) * r_sum_a2'/ sum_a2;
	y0 = (1:K) * c_sum_a2' / sum_a2;
	p0 = angle(sum(sum(z.*a2)));
	c0 = sum(sum(a.*a2)) / sum(a2(:));
	cog = [x0;p0;y0;c0];
else,
	cog = [(M+1)/2;0;(K+1)/2;0];
end;
