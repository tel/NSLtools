global COCHBA;
[L, M] = size(COCHBA);	% p_max = L - 2;

for ch = 1:M,

	p  = real(COCHBA(1, ch));	% order of ARMA filter
	B  = real(COCHBA((0:p)+2, ch));	% moving average coefficients
	A  = imag(COCHBA((0:p)+2, ch));	% autoregressive coefficients
	COCHFIL(:, ch) = freqz(B, A, 2048);
end;
