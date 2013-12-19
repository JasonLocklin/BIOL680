% Week 7 Sandbox

%% Corrcoef()

x = randn(100,4);  % uncorrelated data - 4 signals of 100 points each
x(:,4) = sum(x,2);   % 4th signal is now correlated with first 3
[r,p] = corrcoef(x)  % compute sample correlation and p-values; notice signal 4 is correlated with the others
imagesc(r) % plot the correlation matrix -- note symmetry, diagonal, and r values of ~0.5 for signal 4

%% compute the correlation matrix of a ventral striatal spectrogram
% remember to cd to data first
fname = 'R016-2012-10-03-CSC04a.ncs';
csc = LoadCSC(fname);
 
cscR = Restrict(csc,2700,3300); % risk session only
Fs = 2000;
 
[S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,1:0.25:200,Fs); % spectrogram
 
[r,p] = corrcoef(10*log10(P')); % correlation matrix (across frequencies) of spectrogram
 
% plot
imagesc(F,F,r); 
caxis([-0.1 0.5]); axis xy; colorbar; grid on;
set(gca,'XLim',[0 150],'YLim',[0 150],'FontSize',14,'XTick',0:10:150,'YTick',0:10:150);

% Odd - I get a figure that is much lower in resolution than the figure on the wiki.

% What in this plot changes if you vary the width of the window used to compute the spectrogram? 

% The window size is the implementation of the spectral resolution/noise resiliancy tradeoff

%% Autocorrelation of some synthetic signals
wnoise = randn(1000,1);
[r,p] = corrcoef(wnoise(1:end-1),wnoise(2:end)) % same signal, offset by one same

%% x-corr (correlation accross many lags): White noise
[acf,lags] = xcorr(wnoise,100,'coeff'); % compute correlation coefficients for lags up to 100
plot(lags,acf,'LineWidth',2); grid on;
set(gca,'FontSize',18); xlabel('time lag (samples)'); ylabel('correlation ({\itr})');

%% Now a signal with periodicity
Fs = 500; dt = 1./Fs;
tvec = 0:dt:2-dt;
 
wnoise = randn(size(tvec));
sgnl = sin(2*pi*10*tvec) + wnoise; % 10 Hz sine wave plus Gaussian white noise
 
subplot(211);
plot(tvec,sgnl); grid on;
set(gca,'FontSize',18); title('signal');
 
subplot(212);
[acf,lags] = xcorr(sgnl,100,'coeff');
lags = lags.*dt; % convert samples to time
plot(lags,acf); grid on;
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('autocorrelation');

%% What would you expect the autocorrelation to look like if the above code were modified to use a rectified sine wave (abs(sin(..)))?

wnoise = randn(size(tvec));
sgnl = abs(sin(2*pi*10*tvec)) + wnoise; % 10 Hz sine wave plus Gaussian white noise
 
subplot(211);
plot(tvec,sgnl); grid on;
set(gca,'FontSize',18); title('signal');
 
subplot(212);
[acf,lags] = xcorr(sgnl,100,'coeff');
lags = lags.*dt; % convert samples to time
plot(lags,acf); grid on;
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('autocorrelation');

% A rectifying a sin wave effectively halves it's period and introduces some 
% new components (because it's nolonger a smooth wave). The halved period is
% visible in the autocorrelation, and some roughness is introduced as well. 

% Subtract the mean from the rectified signal and plot the acorr again. Verify that the value at lag 1 matches the output of corrcoef() for lag 1. 

wnoise = randn(size(tvec));
sgnl = abs(sin(2*pi*10*tvec)) + wnoise; % 10 Hz sine wave plus Gaussian white noise
sgnl = sgnl - mean(sgnl)   %Normalize the signal
 
subplot(211);
plot(tvec,sgnl); grid on;
set(gca,'FontSize',18); title('signal');
 
subplot(212);
[acf,lags] = xcorr(sgnl,100,'coeff');
lags = lags.*dt; % convert samples to time
plot(lags,acf); grid on;
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('autocorrelation');


%% Application to real data
% spectrogram with better time resolution (slow!)
Fs = 500;
d = decimate(Data(cscR),4);
[S,F,T,P] = spectrogram(d,hanning(125),115,1:200,Fs);

%% Low gamma autocorrelation
F_idx = find(F > 70 & F < 100); % high-gamma power
F_idx2 = find(F > 50 & F < 65); % low-gamma power
 
pwr = mean(P(F_idx,:)); % average across frequencies
pwr2 = mean(P(F_idx2,:));
 
[ac,lags] = xcorr(pwr2-mean(pwr2),50,'coeff'); % remember to subtract the mean!
lags = lags.*mean(diff(T));
figure;
plot(lags,ac)
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('low-gamma autocorrelation');

% pwr2 is just a bunch of Nan + and Nani's









