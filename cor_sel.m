function [yh, crh] = cor_sel(fname, cr, argin3, sgnw, sw, rw)
% COR_SEL cortical scale-rate weighting
%	yh = cor_sel(fname, cr, NORM, wt);
%	[yh, crh] = cor_sel(fname, cr, [NORM FOUTT FOUTX], sgnw, sw, rw);
%	cor_sel(fname, cr, fnameout, sgnw, sw, rw);
%	fname:		input .cor file
%	cr:			4D cortical representation (sr(up-down)tf)
%	fnameout:	output .cor file
%	NORM:		normalization style, 0<NORM<1, .9 is recommended.
%	FOUTT(X):	overlapped output option, see COR2AUD
%	wt = [sgnw(1)*sw*rw', sgnw(2)*sw*rw'];
%	sgnw:		sign weighting, 2-by-1, ["+" Weight, "-" Weight].
%	sw:			scale weighting, K2-by-1
%	rw:			rate weighting, K1-by-1
%
%	COR_SEL selects rate, scale channels. This function is also good
%	to weight each panel by setting an arbitrary number. If [rw] and
%	[sw] were not given, [sgnw] should be a K2-by-2*K1 matrix. The
%	output can be a file (with FNAMEOUT given) or a variable (with NORM
%	assigned.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 20-Aug-97, non-truncation: FULL, M(4), N(4)
% v1.02: 28-Aug-97, include perfect reconstruction
% v1.03: 04-Sep-97, included FULLOUT.
% v1.04: 03-Oct-97, add causal option
% v1.05: 12-Apr-98, remove non-causal option, SB

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 02-Oct-98, bug: rate filter;
% v1.20: 14-Feb-01, add cr, 4 dimensional representation

% read info 
fcor	= fopen(fname);
[paras, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor);

CR = ~isempty(cr); CROUT = 0;
if (nargout==2 & CR), CROUT = 1; crh = zeros(size(cr)); end;

% rate-scale weighting
if nargin > 4,
	if isempty(sgnw), sgnw = [1; 1];
	else, sgnw = sgnw /sum(sgnw) * 2; end;
	if isempty(sw), sw = sv * 0 + 1;
	else, sw = sw / sum(sw) * K2; end;
	if isempty(rw), rw = rv * 0 + 1;
	else, rw = rw / sum(rw) * K1; end;

	wt = sw(:) * rw(:).';
	wt = [sgnw(1)*wt, sgnw(2)*wt];
else,
	wt = sgnw;
	wt = wt / sum(wt(:)) * (K1 * K2 * 2);
end;

% check weighting matrix size
[k2, k1] = size(wt);
if (k1 ~= K1*2) | (k2 ~= K2),
	fclose(fcor);
	error('Weighting matrix size mismatch !');
end;

% write info
FOUT = (argin3(1)>1);	% note: NORM <= 1

if FOUT,
	fout = fopen(argin3, 'w');
	fwrite(fout, [paras(:); K1; K2; rv(:); sv(:); N(1); M(1); ...
		FULLT; FULLX], 'float');
else,
	STF	= 1000 / paras(1);	% temporal sample frequency
	SRF	= 24;			% spatial sample frequency
	HH	= 0;			% overall response
	Z_cum = 0;			% cumulated reverse response
	NORM = argin3(1);		% normalization style
	if length(argin3) < 2, FOUTT = 0;	else, FOUTT = argin3(2); end;
	if length(argin3) < 3, FOUTX = FOUTT;	else, FOUTX = argin3(3); end;
end;
M(4) = floor(M(1)/2*FULLX);  % dM
N(4) = floor(N(1)/2*FULLT);  % dN

% rate-scale loop
t0 = clock;
for rdx = 1:K1,
	% rate filtering
		if ~FOUT,
			fc_rt = rv(rdx);
			HR = gen_cort(fc_rt, N(2), STF, [rdx K1]);
%			HR = gen_cort(fc_rt, N(2), STF);
		end;
	
	for sgn = [1 -1],

		if ~FOUT,
			% rate filtering modification
			if sgn > 0,
				HR = conj([HR; zeros(N(2), 1)]);
			else,
				HR = [HR(1); conj(flipud(HR(2:N(3))))];
				HR(N(2)+1) = abs(HR(N(2)+2));
			end;
		end;

		for sdx = 1:K2,

			if ~FOUT,
				HS = gen_corf(sv(sdx), M(2), SRF,[sdx K2]);
%				HS = gen_corf(sv(sdx), M(2), SRF);
			end;
	
			% load file
			if CR
				z = squeeze(cr(sdx, rdx+(sgn==1)*K1, :, :));
			else
				z  = corcplxr(fcor, N(1)+2*N(4), M(1)+2*M(4));
			end
			% weighting
			z  = z * wt(sdx, rdx+(1-sgn)/2*K1);
	
			% save file or ifft
			if CROUT,
				crh(sdx, rdx+(sgn==1)*K1, :, :) = z;
			end
			if FOUT,
				corcplxw(z, fout);
			else,
				if NORM >= 0,
					[Z_cum, HH] = corfftc(z, Z_cum, ...
						 N, M, HR, HS, HH);
				else,
					Z_cum = Zcum + z;
				end;
			end;

		end;
	end;
	time_est(rdx, K1, 1, t0);

end;

fclose(fcor);
if FOUT,
	fclose(fout);
	yh = [];
else,
	if NORM >= 0,
		yh = cornorm(Z_cum, HH, N, M, NORM, FOUTT, FOUTX);
	else,
		yh = Z_cum;
	end;
end;

