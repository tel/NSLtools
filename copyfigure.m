function fidnew = copyfigure(varargin)
%
% figID = copyfigure
% figID = copyfigure(figold)
%
% Copies a figure window to a new figure window.

if nargin == 0
	figIDold = gcf;
else
	figIDold = varargin{1};
end

tmpname = 'tempmat';

print(fullfile(tempdir,tmpname), '-dmfile',['-f',num2str(figIDold)]);
addpath(tempdir)
eval(tmpname);
rmpath(tempdir)
delete(fullfile(tempdir,[tmpname,'.m']));
delete(fullfile(tempdir,[tmpname,'.mat']));

if exist('figIDnew')
	figIDnew = gcf;
end
