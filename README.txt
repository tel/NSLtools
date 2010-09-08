%% INTRODUCTION
% the following scripts provide the basic idea of the usage of [nsltools]
% you are welcomed to report any bugs, type-errors as well as suggestions to
% tschi@isr.umd.edu (Taishih Chi), sas@isr.umd.edu (Shihab Shamma)

%% SETUP
% make sure the nsltools is in the path list
% load colormap, filterbank and parameters
loadload;

% you may change the parameters, see WAV2AUD
% paras(1): frame jump, in ms, e.g., 4, 8, 16, 32 ms
% paras(2): Time constant, in ms, e.q., 4, 8, 16, 32, 64 ms
% paras(3): Nonlinear factor, e.g., .1
% paras(4): Octave shift, e.g., -1 for 8k, 0 for 16k (sampling rate)
% rv: rate vector, e.g., 2.^(1:5) or 2.^(1:.25:5); (Hz)
% sv: scale vector, e.g., 2.^(-2:3) or 2.^(-2:.25:3); (cycle/octave)

%% WAVFORM
% sequential format:	see LOADFILE
% .au format:		see AUREAD
% .aiff format:		see AIFFREAD
%x = loadfile('_come.sho');
x = auread('_come.au');	% 8k sample rate

% down sample to 8 kHz if necessary 
% if the sample rate is not the power of 2, you need to use RESAMPLE
% change the paras(4) according to the sample rate
%x = decimate(x, 2); 
%x = x(2:2:length(x));

% if you want to see a schematic plot, you may run
xh = schematc(x);

% make the sequence as unit sequence, i.e., ~ N(0, 1) or
% zero-mean and unit variance (optional but recommended)
% besides, the length 8192, 16384, 32768, etc are preferrable
x = unitseq(x);

% play the sound if you like
soundsc(x, 8000);

%% SPECTRUM
% y is an M-by-128 matrix, M is the # of frames
% a short message will be displayed for every complete octave
y = wav2aud(x, paras);

% plot the spectrogram
aud_plot(y, paras);

%% CORTICAL REPRESENTATION
% static representation (e.g., 60th frame of the spectrogram)
z = aud2cors(y(60, :), sv);
% plot it
cor_plts(z, sv, paras);

% dynamic 4-D representation 
% you can assign an output filename for saving the data see AUD2COR 
% The process will take a while. 15 seconds for 1 second input data on SPARC Ultra-5_10 machine with default parameters
fcorname = '_come.cor';
cr = aud2cor(y, paras, rv, sv, fcorname);

% show the rate-scale energy plot as a function of time
% compute rate-scale-time energy matrix from cortical file
rst = cor2rst(fcorname);

% or computer rate-scale-time energy matrix directly from 4-D variable cr
rst = mean(abs(cr), 4); 

% change back to the jet colormap for displaying purpose 
colormap(jet);
rst_view(rst, rv, sv);


%% MANIPULATIONS 
% at this point, you may do some manipulations
% read Contents.m for Manipulations section

%% SPECTROGRAM RECONSTRUCTION
yh = cor2aud(fcorname);
yh = aud_fix(yh);

%% WAVFORM RECONSTRUCTION
xh = aud2wav(yh, [], [paras 10 1 1]);

