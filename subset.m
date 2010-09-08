function outMatrix = subset(inMatrix, range);
%
%  outMatrix = subset(inMatrix, {range1, range2,...});
%
%    calculates inMatrix(range1,range2,...)
%  
%    useful when inMatrix is the result of a function. The string ':' can
%    be used for any subrange which spans the entire dimension.

outMatrix = inMatrix(range{:});