function num = ext2num(ext)
% EXT2NUM convert SANE Extended to a number
%	num = ext2num(ext);
%	ext: SANE Extended (10 byts)
%	num: floating point number
%	
%	EXT2NUM converts (10-byte) SANE Extendeds to
%	floating point numbers.
% See also: AIFFREAD

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-97

if length(ext) ~= 5, error('Not a Extended number'); end;

% exponent: ext(1);
if ext(1) > 32767,
	ee = ext(1) - 32768
	ss = 1;
else,
	ee = ext(1);
	ss = 0;
end;

% mantissa: ext(2:5);
if ext(2) > 32767,
        fm = ext(2) - 32768;
        ii = 1;
else,
        fm = ext(2);
        ii = 0;
end;

ff = 0;
for k = 5:-1:3,
	ff = (ff + ext(k)) / 65536;
end;

ff = (ff + fm) / 32768;

% number
if ee < 32767,
	num = pow2((-1)^ss*(ii+ff), ee-16383);
elseif ff == 0,
	num = (-1)^ss * Inf;
else,
	numj = NaN;
end;
