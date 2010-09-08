% PD performance measure end
%	pd; 
%
%	PD ends performance measure and reports the various
%	elapsed times.
%	See also: PD, TIME_EST, ALARMSIG

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

% measure the elpased times
%n_flops	= flops - n0_flops;
t_cputime	= cputime - t0_cputime;
t_clock		= etime(clock, t0_clock);

% display
disp(['CPU: ',   num2str(t_cputime), '; ', setstr(9), ...
	'ETIME: ', num2str(t_clock)    ' sec.' ]);
