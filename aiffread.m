function [x, SF] = aiffread(fname);
% AIFFREAD read an AIFF sound file
%	[x, SF] = aiffread(fname);
%	fname	: input AIFF file name
%	x	: audio sample (integer)
%	SF	: sample frequency (float)
%	AIFFREAD reads an AIFF sound file and returns multiple
%	channel audio sample [x] with sample rate [SF].
%	See also: AUREAD, LOADFILE

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 23-Aug-96

% open the AIFF file
fid = fopen(fname);

% FORM file
form_ckid = fread(fid, 4, 'char');		% Form AIFF ckID: 4 char
form_size = fread(fid, 1, 'int32');		% Form AIFF ckSize: int32
form_type = fread(fid, 4, 'char');		% Form AIFF Type: 4 char

% Common Chunk
comm_ckid = fread(fid, 4, 'char');		% Common Chunk ckID: 4 char
comm_size = fread(fid, 1, 'int32');		% Common Chunk ckSize: int32
comm_chno = fread(fid, 1, 'int16');		% Common Chunk numChannel: int16
comm_sfno = fread(fid, 1, 'int32');		% Common Chunk numSampleFrames: int32
comm_sams = fread(fid, 1, 'int16');		% Common Chunk sampleSize: int16
comm_samr = fread(fid, 5, 'uint16');		% Common Chunk sampleRate: extened (10 bytes)
comm_dumm = fread(fid, comm_size-18, 'uint16'); % skip compression fields

% Data Chunk
data_ckid = fread(fid, 4, 'char');		% Sound Data Chunk ckID: 4 char
data_size = fread(fid, 1, 'int32');		% Sound Data Chunk ckSize: int32
data_offs = fread(fid, 1, 'int32');		% Sound Data Chunk offset: int32
data_blks = fread(fid, 1, 'int32');		% Sound Data Chunk blockSize: int32
d_type = ['int' int2str(comm_sams)];
d_size = [comm_chno, comm_sfno];
x = fread(fid, d_size,  d_type);		% Sound Data Chunk blockSize: int8
fclose(fid);

% convert (MAC) extend. to float.
SF = ext2num(comm_samr);