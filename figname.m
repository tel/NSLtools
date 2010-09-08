function figname(NEWNAME) 
% FIGNAME set figure name to machine name or NEWNAME
%	figname;
%	figname(NEWNAME);
%	FIGNAME is handy when you iconify the Matlab figure running
%	on another machine. 

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 12-Aug-97 
% v1.01: 13-Jan-98, changed [unix] to [getenv], shortened machine name

if nargin < 1,
	if isunix,
		R2 = getenv('HOST');
		R1 = findstr('.', R2);
		if length(R1) > 0,
                R2 = R2(1:R1(1)-1);
        end;
		if length(R2) > 4, R2 = R2(1:4); end;
		NEWNAME = [R2 '_' num2str(gcf)];
	else,
		NEWNAME = ['NSL_' num2str(gcf)];
	end;
end;

set(gcf, 'Name', NEWNAME, 'numbertitle', 'off');
