function y = isa1map 
% ISA1MAP true for A1 colormap 
%	y = isa1map;
%	ISA1MAP returns 1 if the current colormap is the A1 colormap and 0 otherwise. 
%	See also: AUD_PLOT, COR_VIEW, AUD2WAV

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 25-Jun-97

% A1 colormap and current colormap
load a1map_a;
currmap = colormap;

% comparison
if size(a1map) == size(currmap),
	if a1map == currmap,
		y = 1;
	else,
		y = 0;
	end;
else,
	y = 0;
end;
