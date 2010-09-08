function savefile(x, fname, prec)
% SAVEFILE save a data file
%	savefile(x, fname, prec);
%	fname: filename with full path and extension
%	prec: precision, e. g., 'short' (default), 'float'

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

fid = fopen(fname, 'w');
if nargin < 3,
	prec = 'short';
end;	
fwrite(fid, x, prec);
fclose(fid);
