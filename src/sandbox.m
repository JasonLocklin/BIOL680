% Sandbox for Week 5 

%% Robust spectral estimation methods

% set up time axis
Fs = 500; % (Hz)
t0 = 0; t1 = 10; % start and end times (s)
tvec = t0:1./Fs:t1; % construct time axis
 
% generate white noise (10 seconds at 500Hz)
x = rand((t1-t0)*Fs,1);
 
% get PSD
wSize = 250; %1/2 second window size seems to work best
nP = 256;
[Porig,Forig] = pwelch(x,rectwin(wSize), wSize/2,nP, Fs);
 
% design filter
W1 = 50; %lower limit (Hz)
W2 = 100; %higher limit (Hz)
nq = 0.5.*Fs;
[b,a] = butter(4,[W1/nq W2/nq]);
% fvtool(b,a); % Plot filter
y = filter(b,a,x);
 
% get PSD
[Pfilt,Ffilt] = pwelch(y, rectwin(wSize), wSize/2,nP, Fs);
 
% plot the resulting PSDs
% subplot(121)
%plot(Ffilt, Pfilt, 'r') % Linear scale
plot(Forig, 10*log10(Porig), 'r'); %Unfiltered in red
hold on;
plot(Ffilt, 10*log10(Pfilt));
grid on;

% Note: press Ctrl+Enter many times to see noise better

%% Ask Matlab for a filter

Wp = [ 50 100] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ 45 105] * 2 / Fs; % stopband
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[b2,a2] = butter(N,Wn); % builds filter
fvtool(b,a,b2,a2) 


%% Greedy

Wp = [ 50 100] * 2 / Fs; 
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[b3,a3] = butter(N,Wn); % builds filter
fvtool(b,a,b3,a3)  %Garbage


%% Greedy with Chebyshev Type I filter
Wp = [ 50 100] * 2 / Fs; 
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b_c1,a_c1] = cheby1(N,0.5,Wn);
fvtool(b2,a2,b_c1,a_c1)   % Note ripple on the passband






