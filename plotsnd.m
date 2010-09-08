function plotsnd(x, SF) 
% PLOTSND plot and play sound 
%	plotsnd(x);
%	plotsnd(x, SF);
%	PLOTSND plots and plays sound

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 14-Aug-97

% plot sound
plot(x);

% interpolation
if nargin < 2,
        SF = 8000;
end;
if SF ~= 8000,
	tic;
	x = interpft(x, round(length(x)*8000/SF));
	%x = resample(x, 8000, SF);
	toc;
end;

% play sound
sound(unitseq(x)/5, 8000);
