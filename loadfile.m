function x = loadfile(fname, prec, n_o_r)
% LOADFILE load a data file
%	x = loadfile(fname);
%	x = loadfile(fname, prec);
%	x = loadfile(fname, prec, n_o_r);
%	fname	: filename with full path and extension
%	prec	: precision, e. g., 'short' (default), 'float'
%	n_o_r	: number of rows. default = 1.
%
%	LOADFILE loads a data file FNAME sequentially in precision 
%	PREC. The data matrix is of the size [n_o_r]-by-??.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

fid = fopen(fname);
if fid ~= -1,
	if nargin < 3,
		n_o_r = 1;
		if nargin < 2,
			prec = 'short';
		end;	
	end;
	x = fread(fid, [n_o_r Inf], prec);
	fclose(fid);
else,
	error('File not found');
end;
