set(0,'DefaultFigureWindowStyle','docked');

%% 1. Compute the PSD of “white noise”, i.e. a random signal 
%  where each sample is drawn independently from the open 
%  interval (0,1) with equal probability. Is it 1/f? How 
%  would you describe its shape? Hint: use the MATLAB function 
%  rand().
rng('shuffle');
d = rand(10000,1); % vector of random numbers [0 1]
Fs = 10;            % fake frequency of fake data
wSize = 1024;        % Window size for 
nP =  1024;         % NFFT size
[Pxx,F] = periodogram(d,rectwin(length(d)), length(d), Fs);
plot(F,10*log10(Pxx));

hold on;
[Pxx,F] =      pwelch(d,rectwin(wSize),wSize/2,nP,Fs);
plot(F,10*log10(Pxx),'r');

% Shape is *not* 1/f, it is completely flat.


%% 2. Compute the PSD of a LFP, simultaneously recorded with 
% the signal above but now from the hippocampus. The ExpKeys 
% specify the filename of a verified hippocampal LFP (in the 
% GoodSWR field, for “good sharp wave-ripple complexes”, 
%     characteristic of the hippocampus). How does it compare 
%     to the ventral striatal LFP?

% V.S. LFP SNIP directly from example %%%%%%%%%%%%%%%%%%%%%%%
clear
cd /home/jalockli/Documents/BIOL680/data/R016-2012-10-08/
csc = LoadCSC('R016-2012-10-08-CSC04d.ncs');
run(FindFile('*keys.m'));
% restrict to prerecord, leaving some time (10s) before rat actually goes on track
csc_pre = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);

csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);
 
% check if sampling is ok
plot(diff(csc_preR)); % only minimal differences
Fs = 1./mean(diff(csc_preR));

% downsample
dsf = 4;
csc_preD = decimate(csc_preD,dsf);
csc_preR = downsample(csc_preR,dsf);
Fs = Fs./dsf;
% Compute the spectrum
wSize = 1024;
[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2, length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             Now, add hippocampus LFP on top:

hold on;
cd /home/jalockli/Documents/BIOL680/data/R016-2012-10-08/
run(FindFile('*keys.m'));
ExpKeys
csc = LoadCSC('R016-2012-10-08-CSC02b.ncs');
csc_preh = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);
csc_preRh = Range(csc_preh);
csc_preDh = Data(csc_preh);

% check if sampling is ok
plot(diff(csc_preRh)); % only minimal differences
Fsh = 1./mean(diff(csc_preRh));
% downsample
dsf = 4;
csc_preDh = decimate(csc_preDh,dsf);
csc_preRh = downsample(csc_preRh,dsf);
Fsh = Fsh./dsf;

wSize = 1024;
%[Pxxh,Fh] = periodogram(csc_preDh,hamming(length(csc_preDh)),length(csc_preDh),Fsh);
[Pxxh,Fh] = pwelch(csc_preDh,hamming(wSize),wSize/2, length(csc_preDh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

% close to 1/f, strong spike at 60Hz, 7Hz activity as well.

%% 3. For both LFPs, explore the effects of the window size 
% parameter of the Welch power spectrum, making sure to vary it 
% across an order of magnitude. What do you notice?
figure;


wSize = 256;
subplot(2,2,1);
[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2, length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on;
[Pxxh,Fh] = pwelch(csc_preDh,hamming(wSize),wSize/2, length(csc_preDh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
title('Window size of 256')

wSize = 512;
subplot(2,2,2);
[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2, length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on;
[Pxxh,Fh] = pwelch(csc_preDh,hamming(wSize),wSize/2, length(csc_preDh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
title('Window size of 512')


wSize = 1024;
subplot(2,2,3);
[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2, length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on;
[Pxxh,Fh] = pwelch(csc_preDh,hamming(wSize),wSize/2, length(csc_preDh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
title('Window size of 1024')


wSize = 2048;
subplot(2,2,4);
[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2, length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on;
[Pxxh,Fh] = pwelch(csc_preDh,hamming(wSize),wSize/2, length(csc_preDh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
title('Window size of 2048')



%% 4. Compare the PSD following the use of decimate() as in the 
% above example, to a PSD obtained from downsampling without 
% using decimate(). Are there any differences? 
figure;

% Decimate (assuming above cells have been run)
subplot(2,1,1);
wSize = 1024;
[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2, length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on;
[Pxxh,Fh] = pwelch(csc_preDh,hamming(wSize),wSize/2, length(csc_preDh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);



%Downsampled (i.e., using csc_preR rather than csc_preD)
subplot(2,1,2);

[Pxx,F] = pwelch(csc_preR,hamming(wSize),wSize/2, length(csc_preR),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
hold on;
[Pxxh,Fh] = pwelch(csc_preRh,hamming(wSize),wSize/2, length(csc_preRh),Fsh);
plot(Fh,10*log10(Pxxh),'r'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

% a mess.
