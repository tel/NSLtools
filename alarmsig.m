function sig = alarmsig
% ALARMSIG alarm signal
%	alarmsig;
%	ALARMSIG plays three-interval-tone signal which can serves 
%	as an alarm 

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 10-Jun-97
% v1.01: 12-Aug-97, replaced [hanning] by expression 

% generating the sound sequence
sig = [sin(2*pi*330*(1:800)/8000) zeros(1, 200) ...
        sin(2*pi*660*(1:800)/8000) zeros(1, 200) ...
        sin(2*pi*440*(1:2000)/8000)];

% windowing and playing the sound
n = 4000;
sig = sig .* (1 - cos(2*pi*(1:n)/(n+1))) / 4;

if nargout == 0,
	sound(sig, 8000);
end;
