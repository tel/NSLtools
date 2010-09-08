% AAA Toolbox -- Acoustic, Auditory, and A1 Cortex.
% Neural Systems Lab., Institutes for Systems Research, UMCP
% Version 1.00 	1-Jun-1997, by Powen Ru (powen@isr.umd.edu)
%	  1.01  18-Aug-1997, added WAV_MOR, SCHEMATC
%
% Acoustic Processing.
%	unitseq	   	  - sequence unification to N(0, 1);
%	aiffread	  - read an AIFF format file
%	ext2num		  - convert SANE Ext. to float.
%	loadfile	  - load a file sequentially
%
% Moving Ripple Generator
%	mvripfft	  - moving ripple generator using FFT
%	mvrip_gui	  - GUI interface for generating moving ripples
%	mvripgui	  - callback function for mvrip_gui interface
%	multimvripfft	  - generalized version of mvripfft
% 
% Auditory Processing.
%	wav2aud		  - auditory spectrogram 
%	wav2y1		  - cochlear filter output
%	aud2wav		  - inverse auditory spectrogram
%	aud2rst	 	  - rate-scale matrix directly from auditory spectrogram
%	aud2sf		  - auditory to scale-frequency energy
%	cochfil		  - information (CF, B, A) of cochlear filters
%	sigmoid		  - nonlinear funcion for cochlear model
%	halfregu	  - half duration regulator 
%	aud_fix		  - fix auditory spectrogrm
%	aud_post	  - auditory post windowing
%	stereaus	  - stereausis
%	
% A1 Cortical Processing.
%	aud2cors	  - static cortical representation
%	cor2auds	  - inverse static cortical representation
%	aud2cor		  - dynamic cortical representation
%	cor2aud		  - inverse dynamic cortical representation
%	cor2rst		  - rate-scale matrix along time axis
%	scalspec	  - scale spectrum
%	gen_cort	  - generate cortical temporal filter
%	gen_corf	  - generate cortical spectral filter
%	corheadr	  - (internal only) read the header of .cor file
%	corcplxr	  - (internal only) read complex numbers from .cor file
%	corcplxw	  - (internal only) write complex numbers into .cor file
%	corfftc		  - (internal only) cortical fft-cumulation routine
%	cornorm		  - (internal only) cortical normalization routine
%
% Manipulations.
%	cor_sel		  - cortical channel selection (weighting)
%	cor_mor		  - cortical morphing
%	cor_map		  - full cortical "mapping"
%	cor_maps	  - static cortical "mapping"
%	aud_maps	  - direct static cortical "mapping"
%	aud_dil		  - dilate auditory spectrogram
%	aud_patch	  - mix timbre-pitch from two auditory spectrograms
%	wav_mor		  - morph two sounds
%
% Utilities.
%	cplxmean	  - re-defined mean for complex number
%	cor_pick	  - pick one channel
%	cor_enhs   	  - enhencement for static cortical representation
%	cor_dist	  - static cortical "distance"
%	cor_info	  - display information of .cor file
%	cochfil		  - cochlear filter bank reader
%	auddist2	  - auditory "distance"
%	rs_sugg		  - rate, scale vectors suggestion
%	rst_view	  - view dynamic cortical representation
%	loadload	  - load necessary things
%	shiftmat	  - shift a matrix	
%
% Graphics / Performance.
%	schematc	  - plot schematic (try this first !)
%	cor_plot	  - plot selected panels
%	cor_plts  	  - plot static cortical representation
%	aud_plot	  - plot auditory spectrogram
%	subplot1	  - custom subplot design
%	cplx_col	  - indices of [a1map] colormap for complex number
%	real_col	  - indices of [a1map] colormap for real number
%	image_c		  - mono-image for complex matrix using arrows
%	image_q		  - mono-image for complex matrix using quivers
%	image_pm	  - mono-image for matrix with positive-negative numbers
%	figname		  - change the name of current figure
%	figsize		  - set current figure size
%	a1fig		  - new figure with A1 colormap and appropriate name
%	isa1map		  - true for the A1 colormap
%	pb, pd		  - computer performance measure
%	time_est	  - estimate elapse time
%	alarmsig	  - alarm signal
 
% Copyright (c) 1997-98 by the Neural Systems Lab., 
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
