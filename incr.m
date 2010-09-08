function incr(incrVar,incrN)
%
% postincr = incr(incrVar) increments incrVar by 1
% postincr = incr(incrVar,incrN) increments incrVar by incrN
%
% see decr

if ~exist('incrN','var'),incrN = 1;end
assignin('base',inputname(1),incrVar+incrN)