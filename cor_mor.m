function yh = cor_mor(fname1, fname2, argin3, par)
% COR_MOR cortical morphing 
%	yh = cor_mor(fname1, fname2, NORM, par);
%	yh = cor_mor(fname1, fname2, [NORM FOUTT, FOUTX], par);
%	cor_mor(fname1, fname2, fnameout, par);
%	fname1, fname2: input files 
%	fnameout: output files
%	NORM: normalization style, see COR2AUD, CORNORM
%	FOUTT(X): overlapped output option, see COR2AUD
%	par: partial of constant morphing, 0 < par < 1; if par < 0, the morphing
%	will take place within the signal.
%	
%	COR_MOR morphs two files according some rules defined in CPLXMEAN.
%	If there is not output file assigned, reconstruction will be executed.
%	The partial PAR is the factor of morphing, e.g., fnameout = 
%	(1-par) * fname1 + par * fname2. The detail method is described
%	in CPLXMEAN. The output can be a file (cortical representation)
%	or a variable (auditory spectrogram).
%	See also: AUD2COR, COR2AUD, CPLXMEAN

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 11-Aug-97, 'float' --> 'int8'
% v1.02: 14-Aug-97, add fclose(fcor(k)), k = 1, 2, instead of fclose(fcor).
% v1.03: 20-Aug-97, add overlapped output (FULLOUT)
% v1.04: 28-Aug-97, include perfect reconstruction
% v1.05: 06-Oct-97, add causal option
% v1.06: 12-Apr-98, remove non-causal option, SB
% v1.07: 15-Apr-98, various length option

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, bugs in various length option 

UNWRAP = 1;
if nargin < 3, argin3 = [1 0 0]; end;
if nargin < 4, par = .5; end;

% read info 
% paras, rv, sv, FULLT, FULLX, should be identical
% length maybe vary
for idx = 1:2,
	fcor(idx) = fopen(eval(['fname' num2str(idx)]));
	[paras, K1, K2, rv, sv, N(:, idx), M, FULLT, FULLX] = corheadr(fcor(idx));
end;

if par >=0 & diff(N(1, 1:2)), INTERP = 1; else, INTERP = 0; end;

if INTERP,
	N(1, 3) = round(par*N(1, 1)+(1-par)*N(1, 2));
	N(2, 3) = 2^nextpow2(N(1, 3)); N(3, 3) = 2 * N(2, 3);
	idx3 = 1:N(1, 3);
	idx1 = linspace(1, N(1, 1), N(1, 3));
	idx2 = linspace(1, N(1, 2), N(1, 3));
else,
	N(:, 3) = N(:, 1);
end;

% write info
FOUT = (argin3(1)>1);	% note: NORM <= 1

if FOUT,
	fout = fopen(argin3, 'w');
	fwrite(fout, [paras(:); K1; K2; rv(:); sv(:); N(1, 3); M(1); ...
		FULLT; FULLX], 'float');
else,
	STF	= 1000 / paras(1); SRF	= 24;
	HH	= 0; Z_cum = 0; NORM = argin3(1);
	if length(argin3) < 2, FOUTT = 0; else, FOUTT = argin3(2); end;
	if length(argin3) < 3, FOUTX = FOUTT; else, FOUTX = argin3(3); end;
end;
M(4) = floor(M(1)/2*FULLX);   % dM
N(4, :) = floor(N(1, :)/2*FULLT);	% dN

% rate-scale loop
t0 = clock;

for rdx = 1:K1,
	% rate filtering
	if ~FOUT,
		fc_rt = rv(rdx);
		HR = gen_cort(fc_rt, N(2, 3), STF, [rdx K1]);
	end;

	for sgn = [1 -1],

		if ~FOUT,
			% rate filtering modification
			if sgn > 0,
				HR = [HR; zeros(N(2, 3), 1)]; HR = conj(HR);
			else,
				HR = [0; conj(flipud(HR(2:N(3, 3))))];
				HR(N(2, 3)+1) = abs(HR(N(2, 3)+2));
			end;
		end;

		for sdx = 1:K2,

			if ~FOUT,
				HS = gen_corf(sv(sdx), M(2), SRF, [sdx K2]); 
			end;

			% load file1 and file2
			z1 = corcplxr(fcor(1), N(1, 1)+2*N(4, 1), M(1)+2*M(4));
			z2 = corcplxr(fcor(2), N(1, 2)+2*N(4, 2), M(1)+2*M(4));

			% morphing 
			if INTERP,
				%z1 = interp1(idx1, z1, idx3);
				%z2 = interp1(idx2, z2, idx3);
				z1 = interp1([1:size(z1,1)], z1, idx1);
				z2 = interp1([1:size(z2,1)], z2, idx2);
			end;
			
			z1 = cplxmean(z1, z2, par, UNWRAP);
	
			% save file or ifft
			if FOUT,
				corcplxw(z1, fout);
			else,
				[Z_cum, HH] = corfftc(z1, Z_cum, N(:, 3), M, HR, HS, HH);
			end;

		end;
	end;
	time_est(rdx, K1, 1, t0);

end;

fclose(fcor(1));
fclose(fcor(2));

if FOUT,
	fclose(fout);
	yh = [];
else,
	yh = cornorm(Z_cum, HH, N(:, 3), M, NORM, FOUTT, FOUTX);
end;
