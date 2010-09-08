function bp_analy(arg1, A)
% BP_ANALY Bandpas filter analysis 
%	bp_analy(h);	% h real, causal
%	bp_analy(H);	% H complex, one-sided
%	bp_analy(B, A);	% B, A cpefficients
SF = 16000;
NF = 8000;
dx = .005;
fontsi = 8;

% calculation
if nargin < 2,	% h or H
	if isreal(arg1),
		h = arg1;
		L = length(h);
		tr = (1:L)/SF*1000;	% ms
		L1 = 2 ^ nextpow2(L);
		L2 = L1 * 2;
		H = fft(h, L2);
		fr = ((1:L1)-1)/L1 * NF;
		absH = abs(H);
        [maxH, maxi] = max(absH);
	else,
		H = arg1;
		L1 = length(H);
		fr = ((1:L1)-1)/L1 * NF;
		h = real(ifft(H, 2*L1));
		absH = abs(H);
	    [maxH, maxi] = max(absH);
		L = round(NF / maxi * 10);
		h = h(1:L) * L1;
		tr = (1:L)/SF*1000;	%ms
	end;
	
	H = H / maxH;
	absH = absH / maxH;
	
	x0 = log2(fr(maxi)/1000);
	xlim = min(1, log2(NF/fr(maxi)));
	xr = (-3:dx:xlim)+x0;
	Hx = interp1(fr(2:L1), 20*log10(absH(2:L1)), 1000*2.^xr);
else,
	B = arg1;
	fr = (1:NF)-1;
	H = freqz(B, A, fr, SF);
	absH = abs(H);
	[maxH, maxi] = max(abs(H));
	H = H / maxH;
	absH = absH / maxH;

	L = round(NF / maxi * 50);
	h = impz(B, A, L) * L; 
	tr = (1:L)/SF*1000;
	xlim = min(1, log2(NF/fr(maxi)));
	x0 = log2(fr(maxi)/1000);
	xr = (-3:dx:xlim)+x0;
	Hx = freqz(B, A, 1000*2.^(xr), SF);
	Hx = 20 * log10(abs(Hx));
end;

% plotting
whitebg('w');
hsub = subplot1(4, 1, .6, .1);

% x - dB
dB3 = -10*log10(2);
xidx = find(abs(diff(Hx > dB3)));
SLOPE = (Hx(xidx+1) - Hx(xidx))/dx;
xedg(1) = interp1(Hx(xidx(1)+[0 1]), xr(xidx(1)+[0 1]), dB3);
xedg(2) = interp1(Hx(xidx(2)+[1 0]), xr(xidx(2)+[1 0]), dB3);
subplot(hsub(1));
plot(xr, Hx, 'r', xedg, [1 1]*dB3, 'b*');
xm = mean(xr);
R1 = 'LR';
for k = 1:2,
	text('position', [x0-1, -20-k*7], ...
		'str', [R1(k) 'S = ' num2str(SLOPE(k)) ' dB/oct'], ...
		'ho', 'le', 'fontwe', 'bold')
end;
set(gca, 'ytick', [-80:20:-20 -3 0], 'fontsi', fontsi);
title(sprintf(['(a) Log. Spectrum: x_0 = %4.2f oct, ', ...
	'BW = %4.2f oct'], x0, diff(xedg)));
R1 = axis;
axis([R1(1:2) -60 5]);
xlabel('Frequency (octave re. 1 kHz)');
ylabel('Log. Magnitude (dB)');
grid on;
drawnow;

% f - amp
fidx = 1:min(maxi*2, length(fr));
BW = 1000 * diff(2.^(xedg));
Q = fr(maxi) / BW;
subplot(hsub(2));
hplot = plot(fr(fidx), absH(fidx), 'g--', fr(fidx), ...
	real(H(fidx)), 'r', fr(fidx), imag(H(fidx)), 'b:');
set(gca, 'ytick', [-.5 0 .5 .7 1], 'fontsi', fontsi, ...
	'yticklabel', ['-.5';' 0 ';'.5 ';'.7 ';' 1 ']);
axis([0 fr(length(fidx))  -1 1]);
h_title = title(...
	sprintf(['(b) Lin. Spectrum: f_0 = %6.1f Hz, BW = %5.1f Hz, ', ...
	'Q = %4.1f'], fr(maxi), BW, Q));
R1 = get(h_title, 'po');
set(h_title, 'po', [R1(1) R1(2)-.1]);
xlabel('Frequency (Hz)');
ylabel('Lin. Amplitude');
grid on;
drawnow;

% f - ph
subplot(hsub(3));
plot(fr(fidx), unwrap(angle(H(fidx)))/pi);
set(gca, 'fontsi', fontsi);
R1 = axis;
axis([0 fr(length(fidx)) R1(3:4)]);
title('(c) Unwrapped Phase Spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase (\pi)');
grid on;
drawnow;

% t - h
subplot(hsub(4));
plot(tr, h);
set(gca, 'fontsi', fontsi);
title('(d) Impulse Response');
xlabel('Time (ms)');
R1 = axis;
axis([0 tr(L) R1(3:4)]);
grid on;

legend(hplot, 'Abs. ', 'Real  ', 'Imag. ', 2);

