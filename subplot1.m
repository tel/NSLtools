function h = subplot1(m, n, sm, sn)
% SUBPLOT1 customizing subplot design
%	h = subplot1(m, n)
%	h = subplot1(m, n, sm, sn)
%	m:  # of plots in a column
%	n:  # of plots in a row
%	sm: vertical space
%	sn: horizontal space
%	h: handle matrix
%
%	SUBPLOT1 generates an array of subplots with flexible space. It
%	returns the (M-by-N) handles of subplot. The unit of SM (SN) is 
%	the height (width) of the axis. If SM and SN are not given, 
%	default .1 will be used.  

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

clf;
h = ones(m, n);

% assume the margin is .2, the space is .1
margin_ratio = 1;

if nargin < 3,
	sm = .1;
	sn = .1;
end;

ht = 1 / (m + 3*margin_ratio*sm + (m-1)*sm);
wd = 1 / (n + 2*margin_ratio*sn + (n-1)*sn);

for mdx = 1:m,
	lo = 1 - (sm*margin_ratio + mdx + (mdx-1)*sm) * ht;
	for ndx = 1:n,
		lf = (sn*margin_ratio + (ndx-1)*(1+sn)) * wd;
		h(mdx, ndx) = subplot('position', [lf, lo, wd, ht]);
	end;
end;

