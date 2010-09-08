function [x,files]= get_timit(spkrid, spchid, files, dr, gd, purp)
% !GET_TIMIT  get time waveform from timit CD on SUN machine
%	[x, files] = get_timit(spkrid, spchid, files, dr, gender, purp)
%	!!! NOT RECOMMENDED FOR EXTERNAL USE
%	x		: speech signal (sampled at 8 KHz),
%	spkrid	: speaker id, 
%	spchid	: speech id,
%	files	: filename of *.wrd files, default = [], (this option protects 
%				klezmer from crashing by using 'unix' command several times.)
%	dr		: dialect region, default = 1,
%	gd		: gender of the speaker, default = 'f',
%	purp	: purpose of the speech, e.g., 'train'(default) or 'test'

% Auther: Taishih Chi(tschi@isr.umd.edu), NSL, UMD
% v1.00: 01-Jun-98

if nargin < 6, purp = 'train'; end,
if nargin < 5, gd = 'f'; end,
if nargin < 4, dr = 1; end,
if nargin < 3, files = []; end,

if isempty(files)
	main_dir = ['/cdrom/cdrom0/timitcd/timit/' purp '/dr'];
	[sta, files] = unix(['ls ' main_dir  num2str(dr) '/' gd '*/*.wrd']);
end

fend = find(abs(files) == 10);
f_no  = length(fend); fend = [0 fend];
fname = files((fend((spkrid-1)*10+spchid)+1):(fend((spkrid-1)*10+spchid+1)-5))

fwrd = fopen([fname '.wrd']);
[ph_ptr, count] = fscanf(fwrd, '%d %d', 2);
ph_sym = fscanf(fwrd, '%s', 1);
spch_bgn = ph_ptr(1);
while count,
	[ph_ptr, count] = fscanf(fwrd, '%d %d', 2);
	if count, spch_end = ph_ptr(2); end 
	ph_sym = fscanf(fwrd, '%s', 1);
end
fclose(fwrd);

fwav = fopen([fname '.wav']);
spch_size = spch_end - spch_bgn;
fseek(fwav, spch_bgn*2+1024, -1); % header length : 1024 bytes
[x, cnt] = fread(fwav, spch_size, 'ushort');
x = rem(x, 256) * 256 + fix(x/256);
x = x - 65536 * sign(sign(x - 32768) +1);
fclose(fwav);

x = decimate(x,2);
