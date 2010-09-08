function D2 = cordist3(z1, z2, SUM, S)
% CORDIST2 Euclidean cortical "distance" (L2)
%	D2 = cordist3(z1, z2);
%	D2 = cordist3(z1, z2, SUM, S);
%
%	z1: complex cortical response matrix 1 (M-by-K)
%	z2: complex cortical response matrix 2 (M-by-K)
%	d:  Euclidean distance
%	SUM: (optional) sum up the distance. 
%
%	See also: AUD2CORS
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 16-Oct-97

if nargin < 3, SUM = 1; end;
if nargin < 4, S = 2; end;

[M, K] = size(z1);


%Kw = max(2.^(-2:.2:3), 1);
%Kw = Kw'/mean(Kw)*K;

D2 = zeros(M, K);

% butter, 2nd, 2/12
%Bx = [0.0495    0.0990    0.0495];
%Ax = [1.0000   -1.2796    0.4776];
% butter, 2nd, 3.12
%Bx = [0.0976    0.1953    0.0976];
%Ax = [1.0000   -0.9428    0.3333];
% butter, 1st, 3.12
%Bx = [0.2929    0.2929];
%Ax = [1.0000   -0.4142];
%%%%% butter, 1st, 2/12
%Bx = [0.2113    0.2113];
%Ax = [1.0000   -0.5774];
Bx = [0.1207   -0.1685    0.2655   -0.1685    0.1207];
Ax = [1.0000   -2.4760    2.8920   -1.6546    0.4185];

ph = (-.8:.05:.8)/2*pi;
cp = cos(ph);
sp = sin(ph);
L = length(ph);
wp = (1+cos(ph-2.3)).^2; wp = wp / sum(wp);

a1 = abs(z1);
a2 = abs(z2);
a_sum = a1 + a2;

% find the COG
[dum, idx] = max(a_sum);
[dum, kc] = max(dum);
mc = idx(kc);

% find local distance
for l = 1:L,

	% real response
	R1 = real(z1)*cp(l) + imag(z1)*sp(l); 
	R2 = real(z2)*cp(l) + imag(z2)*sp(l);

	% HWR
	R1 = max(R1, 0)*2; R2 = max(R2, 0)*2; % better to keep

	% difference
	dR = R1 - R2;

	% balance
	%dR = dR - mean(dR(:)); % better not to have

	%dR = filter(Bx, Ax, dR); % better not to have

	% distance cumulation	
	%D2 = max(D2, (abs(dR).^S)*wp(l));
	D2 = D2 + (abs(dR).^S)*wp(l);	% better
end;

% compute partitioned distance
if SUM*0,
	fprintf(2, '->');
	Mh = fix(M/2);
	Kh = fix(K/2);
	% sum up iso-scale 
	D2p(1) = sum(sum(D2(max(1, mc-Mh):mc, max(1, kc-Kh):kc)));
	D2p(2) = sum(sum(D2(max(1, mc-Mh):mc, kc:min(K, kc+Kh))));
	D2p(3) = sum(sum(D2(mc:min(M, mc+Mh), max(1, kc-Kh):kc)));
	D2p(4) = sum(sum(D2(mc:min(M, mc+Mh), kc:min(K, kc+Kh))));
	D2 = max(D2p);
else,
	fprintf(2, '-x');
	D2 = mean(mean(D2));
end;

D2 = (mean(D2).^(1/S));

