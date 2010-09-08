function decr(decrVar,decrN)
%
% postdecr = decr(decrVar) decrements decrVar by 1
% postdecr = decr(decrVar,decrN) decrements decrVar by decrN
%
% see incr

if ~exist('decrN','var'),decrN = 1;end
assignin('base',inputname(1),decrVar-decrN)