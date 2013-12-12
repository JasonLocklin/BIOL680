
%% plot a simple sinusoid
Fs = 100; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1./Fs:t1; % construct time axis
 
f = 2; % frequency of sine to plot
y = sin(2*pi*f*tvec); % note sin() expects arguments in radians, not degrees (see sind())
 
stem(tvec,y);
