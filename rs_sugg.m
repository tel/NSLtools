function [rv, sv] = rs_sugg(y, paras, SRF, step)
% RS_SUGG rate, scale vector suggestion
%	[rv, sv] = rs_sugg(y, paras, SRF, step);
%	y: auditory spectrum (N-by-M)
%	paras: parameters (see WAV2AUD)
%	step: (optional) step (if scalar) [r_step, s_step] (if vector)  
%
%	RS_SUGG suggests the appropriate rate, scale vector for
%	auditory spectrum Y which was generated according to the 
%	parameter PARAS.
%	See also: WAV2AUD, COR2AUD

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 28-Aug-97

if nargin < 3, SRF = 24; end;
if nargin < 4, step = .5; end;
if length(step) == 1,
	r_step = step;
	s_step = step;
else,
	r_step = step(1);
        s_step = step(2);
end;
 
% dimensions
[N, M] = size(y);

% rate vector
STF = 1000 / paras(1);
rv = 2.^(nextpow2(STF/N)+1:r_step:(nextpow2(STF/2)-1));

% rate vector
%SRF = 24;
sv = 2.^(nextpow2(SRF/M):s_step:(nextpow2(SRF/2)-1));


