function D2 = auddist2(y1, y2)
% AUDDIST2 Euclidean cortical "distance" (L2)
%	D2 = cordist2(y1, y2);
%	y1: complex cortical response matrix
%	y2:  complex matrix
%	d: Euclidean distance
%
%	AUDDIST2 computes the Eulicdean "distance" between two auditory spectrum
%	See also: WAV2AUD
		     
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 16-Oct-97

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, bug: distance for complex numbers 

[N, M] = size(y1);
if M == 128, s1 = std(y1).'; y1 = mean(y1).';
else, s1 = std(y1')'; y1 = mean(y1')'; end;
[N, M] = size(y2);
if M == 128, s2 = std(y2).'; y2 = mean(y2).';
else, s2 = std(y2')'; y2 = mean(y2')'; end;

y1 = y1 - mean(y1);
y2 = y2 - mean(y2);
s  = (s1 + s2) / 2;
dy = y1 - y2;

if min(s),
	dy = dy ./ s;
end;

%D2 = dy.^2;
D2 = dy.*conj(dy);
D2 = sum(D2);
D2 = sqrt(D2);
