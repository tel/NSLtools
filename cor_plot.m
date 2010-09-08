function cor_plot(fname, sgn_sel, r_sel, s_sel, ZLIM, str_desc)
% COR_PLOT cortical plots 
%	cor_plot(fname, sgn_sel, r_sel, s_sel);
%	cor_plot(fname, sgn_sel, r_sel, s_sel, ZLIM)
%	cor_plot(y, para1, rv, sv);
%	cor_plot(y, para1, rv, sv, ZLIM, str_desc);
%
%	USAGE #1:
%	fname:	input .cor file (or y: auditory spectrogram)
%	sgn_sel: sign channel selection, e.g., 1 or 2 ([1 -1]).
%	r_sel:	rate selection, e.g., [3 5 7 9]. (or rv: rate vector)
%	s_sel:	scale selection, e.g., [3 5 7 9 11]. (or sv: scale vector)
%
%	USAGE #2:
%	y: auditory spectrogram
%	para1 = [paras FULLT FULLF]; 
%	paras: see WAV2AUD.
%	%FULLT (FULLX): fullness of temporal (spectral) margin. The value can
%	%	be any real number within [0, 1]. If only one number was
%	%	assigned, FULLT = FULLX will be set to the same value.
%	
%	ZLIM: (optional) magnitude limit. If it is not given, each panel
%		will be normalized by its maximum.
%	str_desc: (optional) string description, This string will appear in
%		the auditory spectrogram plot.
%
%	COR_PLOT generates multiple subplots for a cortical file according
%	various selection. The magnitude (phase) is shown by the saturation 
%	(color). The rate (scale) is in Hz (cycle/octave).
%	IF THE FIRST ARGUEMENT IS AN AUDITORY SPECTROGRAM, this function will
%	plot the decomposed cortical rpresentations directly. Usually, this
%	function is good for quick display. It allows much more freedom to set
%	rate, scale vector. However, it had better to limit the length of
%	vector to be smaller or equal 3 such that all the plots can be
%	accomodated in one landscape page.
%	See also: WAV2AUD, AUD2COR, COR_INFO

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 12-Aug-97, add d_size stuff
% v1.02: 19-Aug-97, include FULL (0 < FULL < 1)
% v1.03: 15-Sep-97, from aud. spec. directly.
% v1.04: 25-Sep-97, added causal option
% v1.05: 12-Apr-98, remove non-causal option

global TICKLABEL;
[N, M] = size(fname);
%load a1map_a;
%colormap(a1map);
ytick = 11:24:128;
if nargin < 5, ZLIM = []; end;
font1 = 6;

%%%%%%%%%%%%%%%%%%%
% PLOT FOR A FILE
%%%%%%%%%%%%%%%%%%%
if N == 1,	% file

	% read info 
	fcor	= fopen(fname);
	[paras, K1, K2, rv, sv, N, M, FULLT, FULLX] = corheadr(fcor);
	L_head  = ftell(fcor);				% length of header

	% dimensions
	dM = floor(M(1)/2*FULLX);
	dN = floor(N(1)/2*FULLT);

	N(2) = N(1) + 2 * dN;
	M(2) = M(1) + 2 * dM;

	if nargin < 2,
		r_sel = 1:K1; s_sel = 1:K2; sgn_sel = 1:2;
	else,
		if max(r_sel) > K1 | max(s_sel) > K2 | max(sgn_sel) > 2,
			error('Size Mismatch !');
		end;
	end;

	Ks1 = length(r_sel);
	Ks2 = length(s_sel);
	Ks3 = length(sgn_sel);

	% log. frequency axis
	y_str = [];
	for fdx = 1:5,
		y_str = [y_str; sprintf('%5.2f', 2^(fdx-3+paras(4)))]; 
	end;

	% allocate subplots
	h = subplot1(Ks1*Ks3, Ks2, .14, .11);

	[h_size, d_size] = corcplxr;
	d_size = N(2)*M(2)*d_size*2 + h_size;
	
	ndx1 = (1:N(1)+2*dN);	% causal case

	for rdxs = 1:Ks1,
		rdx = r_sel(rdxs);
		for sgns = 1:Ks3,
			sgn = sgn_sel(sgns);
			for sdxs = 1:Ks2,
				sdx = s_sel(sdxs);
	
				% seek data panel
				pdx = (((rdx-1) * 2 + sgn-1) * K2 + sdx-1);
				offset = L_head + d_size * pdx; 
				fseek(fcor, offset, -1);
	
				% read file
				z = corcplxr(fcor, N(2), M(2)); 

				% ploting
				if sgn == 1,
					rdx1 = Ks1 - rdxs + 1;
				else,
					rdx1 = rdxs + Ks1*(Ks3-1);
				end;
				subplot(h(rdx1, sdxs));
				if isa1map,
					image(ndx1*paras(1), 1-dM:M(1)+dM, cplx_col(z, ZLIM)');
				else,
					h_img = image_c(z.', [], [32 64], 6, ...
						[[1 length(ndx1)]*paras(1) 1-dM M(1)+dM]);
				end;
				axis xy;

				if FULLT | FULLX, hold on;
					plot([0 N(1)+1 N(1)+1 0 0]*paras(1), ...
				 		[0 0 M(1)+1 M(1)+1 0], 'k--');
				hold off; end;

				% text: rate, scale and maximum magnitude
%				text('position', ...
%					[N(2)*paras(1)/2, ...
%					.9*(M(1)+dM)], ...
%					'str', sprintf(['%5.2f Hz, ', ...
%					'%5.2f c/o (%4.2e)'], ... 
%					rv(rdx)*(3-2*sgn), sv(sdx), ...
%					max(abs(z(:)))), ...
%					'ho', 'ce', 'fontwe', 'bold', ...
%					'fontsi', 10);
				text('position', ...
					[N(2)*paras(1)/2, ...
					.9*(M(1)+dM)], ...
					'str', sprintf(['%5.2f Hz, ', ...
					'%5.2f c/o'], ...
					rv(rdx)*(3-2*sgn), sv(sdx)), ...
					'ho', 'ce', 'fontwe', 'bold', ...
					'fontsi', 10);

				% ylabel, xlabel
				set(gca, 'ytick', ytick, 'fonts', font1);
				if sdxs == round(Ks2/2)+1,
					set(gca, ['y' TICKLABEL], y_str);
				else,
					set(gca, ['y' TICKLABEL], []);
				end;

				if sdxs == 1,
					ylabel('Scale, cyc/oct');
				end;

                if rdxs*(3-2*sgn) ~= 1,
                    set(gca, ['x' TICKLABEL], []);
                end;
                if rdxs*(3-2*sgn) == -Ks1,
                    xlabel('Time (ms)');
                end;
				axis xy; drawnow;
			end;
		end;
	end;
	fclose(fcor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR AN AUDITORY SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else,
	y	= fname;
	M1	= 2^nextpow2(M);	M2 = M1 * 2;
	N1	= 2^nextpow2(N);	N2 = N1 * 2;
	pa	= sgn_sel;
	K1	= length(r_sel);
	K2	= length(s_sel);
	h	= subplot1(2*K1+1, K2+1);
	STF	= 1000/pa(1);

	%if pa FULL

	if nargin < 6, str_desc = []; end;

	% log. frequency axis
	y_str = [];
	for fdx = 1:5,
		y_str = [y_str; sprintf('%5.2f', 2^(fdx-3+pa(4)))];
	end;
	
	%% auditory spec.
	subplot(h(K1+1, 1));
	image((1:N)*pa(1), 1:128, floor(y'/max(y(:))*15.99)+129);
	axis xy;
	
	text('position', [N/2*pa(1), M*.9], 'string', ...
		['Auditory Spectrogram'], ...
		'ho', 'ce', 'fontwe', 'bold', 'fontsi', 10);
	text('position', [N/2*pa(1), N*.1], 'string', str_desc, ...
		'ho', 'ce', 'fontwe', 'bold', 'fontsi', 10);
	set(gca, 'fontsi', 6, 'ytick', ytick, 'linew', 2, ...
		['y' TICKLABEL], y_str, ['x' TICKLABEL], []);
	ylabel('Frequency (kHz)');

	%% pure rate decomposition
	YR = zeros(N1, M);
	for m = 1:M,
		R1 = fft(y(:, m), N2);
		YR(:, m) = R1(1:N1);
	end;
	z = y;  % allocating memory
	for k = 1:K1,
		% decomposition
		fc_rt = r_sel(k);
		HR = gen_cort(fc_rt, N1, STF);

		for m = 1:M,
			R1 = HR .* YR(:, m);
			R1 = ifft(R1, N2);
			z(:, m) = R1(1:N);
		end;

				% plot
		for sgn = [1 -1],
			subplot(h(K1+1-sgn*k, 1));
			if sgn == -1, z = conj(z); end;
			image((1:N)*pa(1), 1:128, cplx_col(z, ZLIM*5)');
			axis xy;
			text('position', [N/2*pa(1), M*.9], 'string', ...
				[num2str(sgn*fc_rt) ' Hz ('  ...
				num2str(round(max(abs(z(:)))*1000)/1000) ...
				')'], 'ho', 'ce', ...
				'fontwe', 'bold', 'fontsi', 8);
			set(gca, 'fontsi', 6, 'ytick', ytick, 'linew', 2, ...
				['y' TICKLABEL], y_str);
			if K1 ~= -sgn*k,
				set(gca, ['x' TICKLABEL], []);
			else,
				xlabel('Time (ms)');
			end;
		end;
		drawnow;
	end;

	%% pure scale decomposition
	clear YR;
	YS = zeros(N, M1);
	for n = 1:N,
		R1 = fft(y(n, :), M2);
		YS(n, :) = R1(1:M1);
	end;

	z = y;	% allocating memory
	for k = 1:K2,
		% decomposition
		fc_sc = s_sel(k);
		HS = gen_corf(fc_sc, M1, 24)';
		for n = 1:N,
			R1 = HS .* YS(n, :);
			R1 = ifft(R1, M2);
			z(n, :) = R1(:, 1:M);
		end;

		% plot
		subplot(h(K1+1, k+1));
		image((1:N)*pa(1), 1:128, cplx_col(z, ZLIM*5)');
		axis xy;
%		text('position', [N/2*pa(1), M*.9], 'string', ...
%			[num2str(fc_sc) ' cyc/oct ('  ...
%			num2str(round(max(abs(z(:)))*1000)/1000) ')'], ...
%			'ho', 'ce', 'fontwe', 'bold', 'fontsi', 8);
		text('position', [N/2*pa(1), M*.9], 'string', ...
			[num2str(fc_sc) ' cyc/oct ('  ...
			num2str(round(max(abs(z(:)))*1000)/1000) ')'], ...
			'ho', 'ce', 'fontwe', 'bold', 'fontsi', 8);

		set(gca, 'fontsi', 6, 'ytick', ytick, 'linew', 2, ...
			['y' TICKLABEL], []);
	end;
	drawnow;

	%% rate-scale decomposition
	Y = zeros(N2, M1);
	for m = 1:M1,
		R1 = fft(YS(:, m), N2);
		Y(:, m) = R1;
	end;
	clear YS;

	z1 = zeros(N, M1);

	for k1 = 1:K1,
		% rate filtering
		fc_rt = r_sel(k1);
		HR = gen_cort(fc_rt, N1, STF);

		for sgn = [1 -1],
			% rate filtering modification
			if sgn > 0,
				HR = [HR; zeros(N1, 1)];
			else,
				HR = [0; conj(flipud(HR(2:N2)))]; 
				%HR(N1+1) = abs(HR(N1+2));
			end;

			for k2 = 1:K2,
				% scale filtering
				fc_sc = s_sel(k2);

				HS = gen_corf(fc_sc, M1, 24);

				% spatiotemporal response
				Z = (HR*HS') .* Y;

				% first inverse fft (w.r.t. time axis)
				for m = 1:M1,
					R1 = ifft(Z(:, m));
					z1(:, m) = R1(1:N);
				end;	% z1: N -by- M1

				% second inverse fft (w.r.t frequency axis)
				for n = 1:N,
					R1 = ifft(z1(n, :), M2);
					z(n, :) = R1(1:M);
				end;	% z: N -by- M

				subplot(h(K1+1-sgn*k1, k2+1));
				image((1:N)*pa(1), 1:128, cplx_col(z, ZLIM)');
				axis xy;
				text('position', [N/2*pa(1), M*.9], ...
					 'string', ...
					[num2str(fc_sc) ' cyc/oct, ' ...
					num2str(fc_rt*sgn) ' Hz. (' ...
					num2str(round(max(abs(z(:)))*1000) ...
					/1000) ')'], 'ho', 'ce', ...
					'fontwe', 'bold', 'fontsi', 8);
				set(gca, 'fontsi', 6, 'ytick', 11:24:128, ...
					['y' TICKLABEL], []);

				if k1*sgn ~= 1,
					set(gca, ['x' TICKLABEL], []);
				end;
				if k1*sgn == K1,
					xlabel('Time (ms)');
				end;
				
				drawnow;

			end;
		end;
	end;
end;
