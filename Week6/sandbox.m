% Week6 tutorial

% help spectrum
%[S,F,T,P] = spectrogram(...) 
%   P is a matrix representing the Power Spectral Density (PSD) of each segment.

%% load the data
% remember to cd to the correct folder here, may need to get this file from the lab database
fname = 'R016-2012-10-03-CSC04a.ncs';
 
csc = LoadCSC(fname);
 
hdr = getHeader(csc);

%crap, I can't run this because I havn't completed Week3's assignment.
