function x = unitseq(x)
% UNITSEQ sequence normalization to N(0, 1)
%	x = unitseq(x);
%	UNITSEQ makes a sequence to be a zero-mean sequence 
%	with variance 1 (i.e., N(0, 1)).

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

x = x - mean(x);	% neutralization
x = x / std(x);		% normalization
