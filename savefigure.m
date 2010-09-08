function savefigure(figNum,figStr)
% SAVEFIGURE(figNum,figNamev)
%
% savefigure saves figure number figNum, using the name figName
% (actually saving 2 files, figName.m and figName.mat). 
%
print('-dmfile',['-f',num2str(figNum)],figStr)

