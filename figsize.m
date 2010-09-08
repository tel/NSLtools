function pixels = figsize(inches)
% FIGSIZE figure size
%	figsize(inches);
%	pixels = figsize(inches);
%	inches: [height width] (inches)
%	pixels: (optional) [height width] (pixels)
%	FIGSIZE sets the size of the figure as desired.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 11-Aug-96

pixels = 97.3 * inches;		% may be machine dependent
pos = get(gcf, 'position');

% chech y dir
screensize = get(0, 'ScreenSize');
if pos(2)+pos(4)+10 > screensize(4),
	pos2 = screensize(4) - pixels(1) - 1;
else,
	pos2 = pos(2)-pixels(1)+pos(4);
end;

set(gcf, 'position', [pos(1) pos2 pixels(2) pixels(1)], ...
	'PaperPositionMode', 'auto');
