function y = halfregu(y)
% HALFREGU half duration regulator 
%	y = halfregu(y);
%
%	HALFREGU is a function which regulates a vector such that
%	the duration is half everywhere
%	See also: WAV2AUD, AUD2WAV
 
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 23-Mar-98
L = length(y);
y = (y > 0);	% y: 0,1 sequence
dy = diff(y);	% up-edge: 1; down-edge: -1.

edge_up = find(dy == 1);
edge_dn = find(dy == -1);

L_dn = length(edge_dn);
L_up = length(edge_up);

if L_dn * L_up,
	dy(edge_dn) = 0*dy(edge_dn);

	if edge_dn(1) > edge_up(1),	% the 1st point is zero
		edge_up(L_up+1) = L;
		edge_dn = round((edge_up(1:L_dn)+edge_up((1:L_dn)+1))/2);
		dy(edge_dn) = dy(edge_dn) - 1;
		y = cumsum([0; dy]);
	else,	% the 1st point is one
		edge_up = [0; edge_up; L];
		edge_dn = round((edge_up(1:L_dn)+edge_up((1:L_dn)+1))/2);
		dy(edge_dn) = dy(edge_dn) - 1;
    	y = cumsum([1; dy]);
	end;
end;
