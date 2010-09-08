function yh = cor_map(fname, argin2, dista, distp, CLIP)
% COR_MAP cortical mapping by cortical "distance" 
%	yh = cor_map(fname, NORM, dista, distp);
%   yh = cor_map(fname, [NORM FOUTT FOUTX], dista, distp, CLIP);
%	cor_map(fname, fnameout, dista, distp);
%   cor_map(fname, fnameout, dista, distp, CLIP);
%	fname	: input .cor file 
%	NORM	: normalization style, see COR2AUD
%	FOUTT(X): overlapped output option, see COR2AUD
%	fnameout: output .cor file
%	dista   : magnitude distance
%	distp	: phase distance
%	NORM	: normalization style, see COR2AUD
%	CLIP	: clipping, see COR_DIST
%	FULLOUT : (optional) overlapped output
%	COR_MAP maps a set of cortical representation into another on the 
%	frame basis. The output can either be afile or a variable.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 20-Aug-97, add non-truncation option, FULL, M(4), N(4)
% v1.02: 28-Aug-97, include perfect reconstruction
% v1.03: 03-Oct-97, add causal option
% v1.04: 12-Apr-98, remove non-causal and SB option

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, bugs: fout, fwrite;

% read info 
fcor	= fopen(fname);
[paras, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor);
if nargin < 5, CLIP = 0; end;		   % default: difference, no clipping
dista = dista';
distp = distp';

% write info
FOUT = (argin2(1)>1) ;	% note: NORM <= 1
if FOUT,
	fout = fopen(argin2, 'w');
	fwrite(fout, [paras(:); K1; K2; rv(:); sv(:); N(1); M(1); ...
		FULLT; FULLX], 'float');
else,
	STF	 = 1000 / paras(1);	  % temporal sample frequency
	SRF	 = 24;				   % spatial sample frequency
	%M(4)= floor(M(1)/2*FULLX);% dM
	%N(4)= floor(N(1)/2*FULLT);% dN
	HH	 = 0;				% overall response
	Z_cum= 0;				% cumulated reverse response
	NORM = argin2(1);		% normalization style
	if length(argin2) < 2, FOUTT = 0;	 else, FOUTT = argin2(2); end;
	if length(argin2) < 3, FOUTX = FOUTT; else, FOUTX = argin2(3); end;
end;
M(4)	= floor(M(1)/2*FULLX);  % dM
N(4)	= floor(N(1)/2*FULLT);  % dN

% rate-scale loop
t0 = clock;
for rdx = 1:K1,
	
	% rate filtering
	if ~FOUT,
		fc_rt = rv(rdx);
		HR = gen_cort(fc_rt, N(2), STF, [rdx K1]);
	end;

	for sgn = [1 -1],
	
		if ~FOUT,
			% rate filtering modification
			if sgn > 0,
				HR = conj([HR; zeros(N(2), 1)]); 
			else,
				HR = [0; conj(flipud(HR(2:N(3))))];
				HR(N(2)+1) = abs(HR(N(2)+2));
			end;
		end;

		for sdx = 1:K2,

			if ~FOUT,
				HS = gen_corf(sv(sdx), M(2), SRF, [sdx K2]);
			end;
	
			% load file
			z  = corcplxr(fcor, N(1)+2*N(4), M(1)+2*M(4));

			% phase shift, magnitude distance
			dp 	= distp(sdx, :);
			da	= dista(sdx, :);

			% transformation
			for n = 1:N,
				% phase
				zp = angle(z(n, :)) + dp;
	
				% magnitude
				za = abs(z(n, :));
				if CLIP,
					za = za .* da;
				else,
					za = za + da;
				end;			 
				z(n, :) = za .* exp(i*zp); 
			end;

			% save file or ifft
			if FOUT,
				corcplxw(z, fout);
			else,
				[Z_cum, HH] = corfftc(z, Z_cum, N, M, HR, HS, HH);
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
	yh = cornorm(Z_cum, HH, N, M, NORM, FOUTT, FOUTX);
end;
