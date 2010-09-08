function y = real_col(y, YLIM)
% REAL_COL indices of A1 colormap for real numbers
%	y = real_col(y);
%	y = real_col(y, YLIM);
%	h = real_col('plot');
%	y	: real matrix --> indices of colormap
%	YLIM	: (optional) level limit
%	REAL_COL converts non-negative real numbers to the indices of 
%	colormap which is used by the cortical representations. If the
%	first input argument is a string, this function will plot the
%	colorbar and return the image handle.
%	See also: CPLX_COL, AUD_PLOT

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 13-Jun-97
% v1.01: 30-Jul-97, make it executable in Matlab4 

% color mapping (while using A1MAP - 256 colors)
%xformed_map = [177 162 147 132 117 102 87 71 ...
%		56 42 27 13 254 238 223 208];
xformed_map = [177 163 151 139 122 105 88 72 ...
                57 42 27 12 253 238 223 208];

global VER;

N_level = length(xformed_map);

if isstr(y),
	y = image(1, 0:1, xformed_map(:)); axis xy;
	if VER < 5,
		set(gca, 'xtick', []);
	else,
		set(gca, 'xtick', [], 'yaxisloc', 'right');
	end;
	ylabel('Level (normalized)');
else,
	% clipping
	if nargin < 2, YLIM = max(y(:));
	elseif isempty(YLIM), YLIM = max(y(:)); end;

	% negative value ?
	if min(y(:)) < 0,
		y = min(max(y, -YLIM), YLIM);
		y = y - min(y(:));
		y = y / max(y(:)) * YLIM;
	end;

	% color coding
	y = min(max(ceil(real(y)/YLIM*N_level),1),N_level);
	if VER < 5,
		for n = 1:size(y, 1),
			y(n, :) = xformed_map(y(n, :));
		end;
	else,
		y = xformed_map(y);
	end;
end;
