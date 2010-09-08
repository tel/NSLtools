function h = a1fig(arg);
% A1FIG new figure with A1 colormap and appropriate name
%	h = a1fig;
%	a1fig(fig_handle);
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 30-Oct-97 
% v1.01: 13-Jan-98, changed [unix] to [getenv], shortened machine name

global a1map;

if isunix,
	%[R1, R2] = unix('echo $HOSTNAME');
	R2 = getenv('HOST');
	R1 = findstr('.', R2);
	if length(R1) > 0,
		R2 = R2(1:R1(1)-1);
	end;
	if length(R2) > 4,
		R2 = R2(1:4);
	end;
else,
	R2 = 'NSL';
end;

if nargin,
	set(arg, 'Name', [R2 '_' num2str(arg)], 'Colormap', a1map, ...
		'numbertitle', 'off');
else,
	h = figure('Colormap', a1map, 'numbertitle', 'off');
	set(h, 'Name', [R2 '_' num2str(h)]);
end;
