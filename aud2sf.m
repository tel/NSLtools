function z = aud2sf(y, sv);
% AUD2SF auditory to scale-frequency energy
%	z = aud2sf(y, sv);
%	y   : auditory spectrogram, N-by-M, where
%	sv  : scale vector in cyc/oct, e.g., 2.^(-2:.5:3).
%
% AUD2SF 
% See also: WAV2AUD.

z = 0;
N = size(y, 1);

t0 = clock;
for n = 1:N,
	z = z + abs(aud2cors(y(n, :), sv));
	time_est(n, N, 50, t0);
end;

z = z / N;
