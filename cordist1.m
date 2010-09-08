function D2 = cordist1(z1, z2, SUM, S)
% CORDIST2 Euclidean cortical "distance" (L2)
%	D2 = cordist2(z1, z2);
%	D2 = cordist2(z1, z2, SUM, S);
%
%	z1: complex cortical response matrix 1 (M-by-K)
%	z2: complex cortical response matrix 2 (M-by-K)
%	d:  Euclidean distance
%	SUM: (optional) sum up the distance. 
%
%	CORDIST2 computes the Eulicdean "distance" between two static cortical
%	representations.
%	See also: AUD2CORS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 16-Oct-97

if nargin < 3, SUM = 1; end;
if nargin < 4, S = 2; end;

[M, K] = size(z1);

a1 = abs(z1);
a2 = abs(z2);

[p1, x1] = max(a1);
[p1, y1] = max(p1);	% best scale
[a1, x1] = max(a1(:, y1)); % best response, best frequency
p1 = angle(z1(x1, y1)); % best phase;

[p2, x2] = max(a2);
[p2, y2] = max(p2); % best scale
[a2, x2] = max(a2(:, y2)); % best response, best frequency
p2 = angle(z2(x2, y2)); % best phase;

D2 = (abs((x1-x2))^S+...
	abs((y1-y2))^S+...
	abs((angle(exp(i*(p1-p2)))))^S+...
	abs((a1-a2)/1000000).^S).^(1/S);
