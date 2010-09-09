function nslmovie(rst, rv, sv, paras, sound, samprate, outfile)
%% NSLMOVIE avi movie of rate-scale plot with sound
%    rst      : scale-rate-time matrix
%    rv       : rate vector in Hz
%    sv       : scale vector in cyc/oct
%    paras    : NSLTOOLS parameter vector
%    sound    : time representation of analyzed sound
%    samprate : sampling rate for `sound`
%    outfile  : filename for output of avi movie
  
% Author: Joseph Abrahamson (tel@jhu.edu)
% Depends: Malcolm Slaney, "QuickTime tools for Matlab," Interval Technical Report 
%          1999-066, 1999. <http://cobweb.ecn.purdue.edu/~malcolm/interval/1999-066/>

%% Code
    
    [ns, nr, nframe] = size(rst);
    nr = nr/2;
    if nr ~= length(rv)
        error('Scale key vector does not correspond to RST matrix');
    end
    
    if ns ~= length(sv)
        error('Rate key vector does not correspond to RST matrix');
    end

    p = struct;                         % Function parameter
                                        % structure. NOTE that it is not
                                        % the same as `paras`.
    p.blur = 5;                         % Smoothing factor

    p.slowdown = 2^0;                     % Slowdown factor on the replay
    p.block    = 200;                   % Number of samples to block
                                        % during slowdown.
    
    
    
    %% Initialize movie    
    MakeQTMovie('start', outfile);
    MakeQTMovie('quality', 1);
    MakeQTMovie('size', [640, 320]);
    
    
    %% Determine framerate
    
% $$$     timelen = length(sound)/samprate    % length in sec
% $$$     nframes = timelen ...
% $$$               * 1000 ...                % 1000 ms/s
% $$$               / paras(1);               % ms/frame
% $$$     
% $$$     framerate = nframes/timelen;
    
    % I HAVE NO IDEA WHY THIS IS WORKING: samprate/8000!!!
    MakeQTMovie('framerate', ceil(samprate/8000*1000/paras(1)/p.slowdown));
    
    
    %% Build movie frames

    
    % Each frame, rst(:, :, i), needs to be transformed to be
    % viewable similar to what rst_plot produces. The final axes will be 
    
    % y+ : scale increasing
    % x+ : rate increasing (starts with negative mirror)
    % t+ : time increasing
    
    % but the RST matrix itself has a structure more like
    
    % (scale increasing)
    % x [(rate decreasing below 0) + (rate increasing above 0)] 
    % x (time)
    
    % So to get imagesc to produce the proper orientation the RST matrix
    % must be split into two halves along the scale dimension, the left
    % half needs to be transposed, the right half needs to be flipped
    % vertically, then they are re-adheared. 
    
    % Additionally, each half undergoes a separate cubic interpolation to
    % increase the resolution by the parameter `p.blur`.
    
    % Axis scaling constants
    maxr = max(rv); stepr = maxr/length(rv);
    rlabs = [-rv(end:-1:1) 0 rv(1:end)];
    maxs = max(sv); steps = maxs/length(sv);
    slabs = sv(end:-1:1);
    
    % Make frame rescale function
    rstmax = max(max(max(rst)));
    rstmin = min(min(min(rst)));
    rescale = @(frame) ((frame-rstmin)./rstmax)*5;

    percentstep = 2;
    
    for i = 1:nframe

        %% Generate each image frame
        
        % Split halves while reversing
        left = rst(ns:-1:1, nr:-1:1, i);
        right = rst(ns:-1:1, (nr + 1):2*nr, i);

        % Smooth each half
        left = interp2(left, p.blur, 'cubic');
        right = interp2(right, p.blur, 'cubic');
        
        [ysz, xsz] = size(left);
        xsz = xsz*2+1;                  % Includes the 0 column
        
        % Adhear the two halves including blank center column
        frame = [rescale(left) zeros(ysz, 1), rescale(right)];
        
        image(frame);
        
% $$$         % Fix the axes labels
% $$$         xstep = xsz/length(rlabs);
% $$$         % The 1/2 constant shifts the ticks to the center of pixels
% $$$         set(gca, 'XTick', 1/2+xstep/2:xstep:xsz);
% $$$         set(gca, 'XTickLabel', rlabs);
% $$$ 
% $$$         ystep = ysz/length(slabs);
% $$$         set(gca, 'YTick', ystep/2:ystep:ysz);
% $$$         set(gca, 'YTickLabel', slabs);
        
        %% Save the image frame to the movie
        
% $$$         MakeQTMovie('addfigure');
        MakeQTMovie('addmatrix', frame);
    
        
        %% Output some timing info

        if i/nframe*100 >= percentstep
            fprintf('(%d%%) ', percentstep)
            percentstep = percentstep + 2;
        end
        
    end
    
    %% Finalize movie

    % Pitch-constant sound stretching
    if p.slowdown > 1
        so = [];
        for i = 0:length(sound)/p.block-1
            for j = 1:p.slowdown
                so = [so, sound(1+(i*p.block):(i+1)*p.block)'];
            end
        end
    else
        so = sound;
    end
    
    % Add the sound track to the movie
    MakeQTMovie('addsound', so, samprate);
    
    MakeQTMovie('finish');
    MakeQTMovie('cleanup');
