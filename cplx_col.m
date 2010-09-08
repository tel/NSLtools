function y = cplx_col(z1, z2, N_level, N_color, ZLIM)
% CPLX_COL complex number display preprocessor.
%	c = cplx_col(z); (default)
%	c = cplx_col(z, ZLIM);
%	c = cplx_col(mag, pha ,N_level, N_color);
%	c = cplx_col(mag, pha ,N_level, N_color, ZLIM);
%	h = cplx_col('plot');
%	h = cplx_col('plot', N_level, N_color);
%	z = mag * exp(i * pha): complex matrix
%	N_level: # of magnitude levels
%	N_color: # of colors
%	ZLIM: magnitude limit 
%
%	CPLX_COL displays a complex matrix by using a special color map 
%	that has N_LEVEL levels of saturation for N_COLOR colors. if the
%	first argument is 'plot', it will plot the colormap.

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97
% v1.01: 14-Jun-97, include colorbar function
% v1.02: 18-Aug-97, include ZLIM < 0 (autoscaling) option.
% v1.03: 28-Oct-97, remove bug: when ZLIM == 0
% v1.04: 11-Nov-97, bug in h = cplx_col('plot', N_level, N_color);

if isstr(z1),
	if nargin < 3, N_level = 16; N_color = 16;
	else, N_color = N_level;  N_level = z2; end;
	R1 = 1:N_level*N_color;
	y = image(-1:1, 0:1, reshape(R1, N_level, N_color)); axis xy;
	xlabel('Phase (\pi)'); ylabel('Level (normalized)');
	title('Complex Number Mapping');
else,

	% phase -pi < pha < pi
	% magnitude 0 < mag
	if nargin < 3,
		mag = abs(z1);
		pha = angle(z1);
		if nargin == 2,
			if isempty(z2),
				ZLIM = max(mag(:));
			elseif z2 < 0,
				ZLIM = max(mag(:));
			else,
				ZLIM = z2;
			end;
		else,
			ZLIM = max(mag(:));
		end;
		N_level = 16;
		N_color = 16;

	else,
		mag = z1;
		pha = z2;
		if nargin < 5,
	   		ZLIM = max(mag(:));
		elseif isempty(ZLIM),
			ZLIM = max(mag(:));
		elseif ZLIM < 0,
			ZLIM = max(mag(:));
	   	end;
	end;

	% magnitude coding, 0->ZLIM ===> 1:N_level
	if ZLIM == 0, ZLIM = 1; end;
	mag = mag / ZLIM;			   		% 0->1 
	mag = max(min(ceil(mag*N_level),N_level),1);	% 1:N_level

	% phase coding, -pi->pi ===> 0:N_color-1
	pha = (pha / pi + 1) / 2;			% 0->1
	pha = max(min(ceil(pha*N_color),N_color),1)-1;	% 0:N_color-1

	% index matrix
	y =  pha * N_level + mag;

end;
