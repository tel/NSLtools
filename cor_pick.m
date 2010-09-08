function z = cor_pick(fname, sgn_sel, r_sel, s_sel)
% COR_PICK cortical plots 
%	cor_pick(fname, sgn_sel, r_sel, s_sel);
%	fname:	 input .cor file
%	sgn_sel: sign selection, e.g., 1 or 2 (for 1, -1).
%	r_sel:	 rate selection, --> rv(r_sel)  
%	s_sel:	 scale selection, --> sv(s_sel) 
%		will be normalized by its maximum.
%
%	COR_PICK picks the cortical representation for [r_sel]-th rate and
%	[s_sel]-th scale panel. Note: only one panel is allowed.
%	See also: AUD2COR

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 12-Aug-97
% v1.01: 19-Aug-97, add FULL panel option
% v1.02: 12-Apr-98, remove non-causal option

% read info 
fcor    = fopen(fname);
[paras, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor);
L_head  = ftell(fcor);				% length of header

% dimensions
N = N(1) + 2 * floor(N(1)/2*FULLT);
M = M(1) + 2 * floor(M(1)/2*FULLX);

if r_sel > K1 | s_sel > K2,
	error('SelectionExceedLimit !');
end;

% seek data panel
pdx = (((r_sel-1) * 2 + sgn_sel-1) * K2 + s_sel-1);
[h_size, d_size] = corcplxr;
offset = L_head + (N * M * d_size * 2 + h_size) * pdx; 
fseek(fcor, offset, -1);

% read file
z  = corcplxr(fcor, N, M); 
fclose(fcor);
