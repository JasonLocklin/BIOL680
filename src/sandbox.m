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

%% Phase responses and filtfilt()


Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
s1 = sin(2*pi*80*tvec+pi/6);
s2 = sin(2*pi*40*tvec);
s = s1 + s2;
 
sf = filter(b_c1,a_c1,s);
 
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);


%% Run fvtool again on the Butterworth and Chebyshev filters above, 
% and now select the Phase Response button in the top left of the window. 
fvtool(b2,a2,b_c1,a_c1);
legend({'Butterworth', 'Chebyshev'})
% Click "Analysis" -> "Phase Response"

% They are simmilar for Normalized frequency < 0.4 rad/sample, then
% diverge. Butterworth is always better.

%% Filter signal forward, then backward, creating an overall null shift
sf = filtfilt(b_c1,a_c1,s);
 
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);

%% compare freq responses
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
x = rand(size(tvec)); % white noise input
[P,F] = pwelch(x,hanning(512),256,2^14,Fs);
 
y1 = filter(b_c1,a_c1,x);
[P1,F1] = pwelch(y1,hanning(512),256,2^14,Fs);
 
y2 = filtfilt(b_c1,a_c1,x);
[P2,F2] = pwelch(y2,hanning(512),256,2^14,Fs);
 
plot(F,10*log10(P),F,10*log10(P1),F,10*log10(P2));
legend({'original','filter','filtfilt'});


%% Some typical neuroscience applications
%% Removing 60Hz line noise in the data
[b,a] = butter(10, [59 61] * 2 / Fs, 'stop');
fvtool(b,a);

%Better?
Wp = [ 58 62] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ 59 61] * 2 / Fs; % stopband
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[b2,a2] = butter(N,Wn, 'stop'); % builds filter
fvtool(b2,a2) 


%% Another method (best)
[z,p,k] = butter(10, [59 61] * 2 / Fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
h = dfilt.df2sos(sos,g); % create filter object
fvtool(h);

% test with some white noise
x_pre = rand((t1-t0)*Fs +1,1);
tvec = t0:1./Fs:t1; % construct time axis

hum = 0.5* sin(2*pi*60*tvec);
x = hum(:)+x_pre % Add some 60hz hum to the signal
wSize = 256; %1/2 second window size seems to work best
nP = 256;
[Porig,Forig] = pwelch(x,rectwin(wSize), wSize/2,nP, Fs);

%Test this notch filter:
sf = filtfilt(sos, g,x);
% get PSD
[Pfilt,Ffilt] = pwelch(sf, rectwin(wSize), wSize/2,nP, Fs);
 
% plot the resulting PSDs
% subplot(121)
%plot(Ffilt, Pfilt, 'r') % Linear scale
plot(Forig, 10*log10(Porig), 'r'); %Unfiltered in red
hold on;
plot(Ffilt, 10*log10(Pfilt));
grid on;



%% Detecting movement artifacts


cd /home/jalockli/Documents/BIOL680/data/R016-2012-10-08/
csc = LoadCSC('R016-2012-10-08-CSC02b.ncs');
cscR = Restrict(csc,1270,1272);
plot(Data(cscR))


%% Guess the frequency
x = Data(cscR);
tvec = Range(cscR);
 
Fs = 2000;
Wp = [ 180 220] * 2 / Fs;
Ws = [ 178 222] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); % determine filter parameters
[b_c1,a_c1] = cheby1(N,0.5,Wn); % builds filter
 
%fvtool(b_c1,a_c1); % remember to check your filter!
 
y = filtfilt(b_c1,a_c1,x);
plot(tvec,x,'b',tvec,y,'r');


chew_power = y.^2;


chew_power_filtered = medfilt1(chew_power,101); % filter window is specified in samples, so this is ~50ms
plot(tvec,x,'b',tvec,chew_power_filtered,'r');



%% Assignment
cd /home/jalockli/Documents/BIOL680/data/R042-2013-08-18/



% extract varargins
 
% load signal
csc = LoadCSC('R042-2013-08-18-CSC03a.ncs');
cscR = Restrict(csc,4*10^3, 4.01*10^3);
plot(Data(cscR))


% filter in the ripple band
 
% convert to power envelope
 
% convert to z-score (SDs from the mean; use nanmean() and nanstd())
 
% find times when above threshold
 
% find crossings from below to above threshold and vice versa (can use diff())
 
% get center time and power, return
 








































