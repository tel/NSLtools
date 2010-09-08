function swapbyte(fnamein, fnameout)
% SWAPBYTE to swap high/low byte
%	swapbyte(fname);
%	swapbyte(fnamein, fnameout);
%	SWAPBYTE swaps high/low byte for [fnamein]. The new file will be
%	[fnameout]. If the output filename is not given, the input file will
%	be overwriten. This utiliy is good to convert files from other OS.
%	See also: LOADFILE

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 20-Oct-97

if nargin < 2, fnameout = fnamein; end;

fid = fopen(fnamein);
if fid ~= -1,
        [x, count] = fread(fid, Inf, 'int8');
        fclose(fid);

	if rem(count, 2),
		error('File size not an even number');
	else,
		x = reshape(x(:), 2, count/2);
		x = x([2 1], :);
		fid = fopen(fnameout, 'w');
		fwrite(fid, x, 'int8');
		fclose(fid);
	end;
else,
        error('File not found');
end;

