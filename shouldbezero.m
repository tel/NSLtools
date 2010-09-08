function out = shouldbezero(test, epsilon)
% SHOULDBEZERO complains when given a number that should be small but isn't.
%
%       shouldagain = shouldbezero(test, epsilon)
%
%       SHOULDBEZERO is an M-file that takes test should to epsilon
%       and beeps and prints an error message. It is most useful for
%       debugging, especially for testing the imaginary part of an
%       inverse fourier transform that should be completely real but
%       usually isn't due to numerical imprecision.
%

out = test;
if abs(test) > abs (epsilon)
 fprintf(1, ...
  '%c\n%g should be zero, but its magnitude is greater than %g!\n\n' ...
  ,7, test, epsilon);     % 7 = ^G = 'beep'
end


