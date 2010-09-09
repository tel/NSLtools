function tomov(wavfile, f0, outfile)
%% TOMOV create avi movie of long clip of a song
%  wavfile : path to wav file for conversion
%       f0 : sampling frequency of input wav file
%  outfile : path to output mpg file

% TOMOV is an automated script for the generation of movies of the
% cortical response to an input sound. The bulk of the work is similar to
% the smaller, batch function NSLMOVIE, but the song is processed in a
% number of short chunks in order to maintain a lower memory
% overhead. Additionally, the entire auditory-cortical conversion is run
% automatically with a default set of parameters.
  
% Author: Joseph Abrahamson (tel@jhu.edu)
% Depends: Malcolm Slaney, "QuickTime tools for Matlab," Interval Technical Report 
%          1999-066, 1999. <http://cobweb.ecn.purdue.edu/~malcolm/interval/1999-066/>

loadload;

%% TOMOV internal parameters

slowdown = 1;                           % Framerate adjustment. Hopefully
                                        % not necessary.
nblur = 2;                              % Smoothing interpolation scale

%% Set up some basic NSLTOOLS parameters

paras = [8, ...                         % 8ms frame skips
         8, ...                         % 8ms time constant
         -1, ...                        % Linear compression at first
         0];                            % 16k sampling rate.

f1s = 16000;                             % To match the paras sampling rate

% Rate and scale vectors. Sort of arbitrary right now.
sv = (1.3).^(-4:11); nsv = length(sv);
rv = (1.4).^(-3:8); nrv = length(rv);

% Labels and ticks for later (UNUSED right now)
maxr = max(rv); stepr = maxr/length(rv);
rlabs = [-rv(end:-1:1) 0 rv(1:end)];
maxs = max(sv); steps = maxs/length(sv);
slabs = sv(end:-1:1);


%% Adjust the wavfile for processing

% Read wav file, monaural
x = wavread(wavfile);
x = sum(x, 2); 

% Adjust to be ~ N(0, 1)
x = unitseq(x);

% Downsample to 16k, divide by 3000 for a simpler fraction
p  = ceil(f1s/3000);
q  = ceil(f0/3000);
x  = resample(x, p, q);
f1 = ceil(f0 * p/q);

% Some global measures of the sound
xl = length(x);


%% Basic processing parameters

chunksize = 3 * f1; % samples of song to process at a time
nchunks = ceil(xl/chunksize);

movwidth = 640; %px

MakeQTMovie('start', outfile);
MakeQTMovie('quality', 1);

% equal aspect, rv:sv :: w:h
MakeQTMovie('size', [movwidth, ceil(movwidth*(nsv/nrv/2))])

% framerate derived in NSLMOVIE.m
MakeQTMovie('framerate', ceil(f1/f1s * 1000/paras(1))/slowdown);


%% Chunk processing

for ci = 1:nchunks
    
    % Get the chunk of the wav file
    if ci ~= nchunks
        xc = x(((ci-1)*chunksize+1) : (ci*chunksize+1));
    else
        % it's the last chunk so the second index is `end`
        xc = x(((ci-1)*chunksize+1) : end);
    end

    % Output timing diagnostics
    fprintf('Processing chunk: %d of %d... \n', ci, nchunks);
    
    % Get the auditory representation
    % (includes cube-root compression)
    y = power(wav2aud(xc, paras), 1/3);
    
    % Compute the RST-matrix without using a temp file
    cr = aud2cor(y, paras, rv, sv, 'tmpxxx');
    clear y;
    rst = mean(abs(cr), 4);
    clear cr;
    
    % Get tick sizes, #rates #scales
    [ns, nr, nframe] = size(rst);
    nr = nr/2;
    
    rstmax = max(max(max(rst)));
    rstmin = min(min(min(rst)));
    rescale = @(frame) ((frame-rstmin)./rstmax)*5;
    

    for fi = 1:nframe
        
        % Separate halves and smooth
        left = rst(ns:-1:1, nr:-1:1, fi);
        right = rst(ns:-1:1, (nr + 1):2*nr, fi);
        left = interp2(left, nblur, 'cubic');
        right = interp2(right, nblur, 'cubic');
        
        [ysz, xsz] = size(left);
        xsz = xsz*2+1;                  % Includes the 0 column
        
        % Adhear the two halves including blank center column
        frame = [rescale(left) zeros(ysz, 1), rescale(right)];
        
        % Write it to the movie
        MakeQTMovie('addmatrix', frame);
    
    end
    clear rst;
   
end


%% Movie finalization

% Add the sound track to the movie
x = resample(x, 1, slowdown);
MakeQTMovie('addsound', x, f1/slowdown);

MakeQTMovie('finish');
MakeQTMovie('cleanup');


