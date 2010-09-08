function time_est(n, N, LAP, t0, COMMENT)
% TIME_EST estimate elapsed time 
%	time_est(n, N, LAP, t0);
%	time_est(n, N, LAP, t0, COMMENT);	
%	n	: current loops
%	N 	: total number of loops
%	LAP	: the lap to display the message
%   t0	: phase distance
%	COMMENT: (optional) any string
%
%	TIME_EST estimate elapse time and display the message
%	every LAP loops. It is good for users to know how the
%	program progressed.
%	See also: PB, PD

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

if ~rem(n, LAP),
	t(1) = etime(clock, t0);
	t(2) = t(1) * N /n;
	t(3) = t(2) - t(1);
	s = rem(t, 60);
	m = floor(t/60);
	h = floor(m/60);
	m = rem(m, 60);
	if nargin < 5, COMMENT = ''; end;

	fprintf(['%5.1f%% processed; elaps. %d:%d:%4.1f; ' ...
		'esti. %d:%d:%4.1f; rem. %d:%d:%4.1f; %s\n'], ...
		n/N*100, h(1),m(1),s(1), h(2),m(2),s(2), h(3),m(3),s(3), COMMENT);	
end;
