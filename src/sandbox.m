set(0,'DefaultFigureWindowStyle','docked');

% first, cd to where the data you just grabbed is located
cd /home/jalockli/Documents/BIOL680/data/R042-2013-08-18/

%% load the data (note, may need to unzip position data first)
fc = FindFiles('*.t');
S = LoadSpikes(fc);

[csc,csc_info] = LoadCSC('R042-2013-08-18-CSC03a.ncs');

[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );


%% Week 2 assignment

% Creating a 100Hz signal
fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz
 
freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal
 
ax1 = subplot(211);
stem(tvec,y); title('original'); %Stem and circle plot

%Downsample (every fourth sample)
subsample_factor = 10; %4 Produces a reasonable result, 10 results in garbage
 
tvec2 = tvec(1:subsample_factor:end); % take every 4th sample
y2 = y(1:subsample_factor:end);
 
ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled'); 
xlabel('time (s)');

% Rescale axis
xl = [1 1.04];
set(ax1,'XLim',xl); set(ax2,'XLim',xl);

% Recover origional signal
hold on;
 
y_interp = interp1(tvec2,y2,tvec,'linear');
p1 = plot(tvec,y_interp,'b');
 
y_interp2 = interp1(tvec2,y2,tvec,'nearest');
p2 = plot(tvec,y_interp2,'g');

% Not able to plot spline because of this error
% Error using griddedInterpolant/subsref
% LAPACK loading error:
% dlopen: cannot load any more object with static
% TLS
% 
% Error in interp1 (line 192)
%         VqLite = F(Xqcol);
 
legend([p1 p2],{'linear','nearest'},'Location','Northeast'); legend boxoff

%
sin(2*pi) == 0 % ...right? The answer might surprise you.
fprintf('Welcome to numerical computing!\n');

%% Just at the Nyquist frequency
figure(2)
subsample_factor = 10;
 
tvec2 = tvec(2:subsample_factor:end); % best case scenario -- can detect 100Hz signal
y2 = y(2:subsample_factor:end);
 
subplot(212)
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% 2kHz Fs, 100Hz signal with 450Hz signal superimposed
figure(3)
 
fs1 = 2000;
tvec = 0:1/fs1:4;
 
freq1 = 100;
freq2 = 450; % note, above Nyquist frequency for our target subsampled Fs
 
y = sin(2*pi*freq1*tvec) + 0.5.*sin(2*pi*freq2*tvec);
 
subplot(211)
stem(tvec,y)
set(gca,'XLim',xl);

%% ss -- we don't care about the 450Hz signal, but...
subsample_factor = 4;
 
tvec2 = tvec(1:subsample_factor:end);
y2 = y(1:subsample_factor:end);
 
subplot(212)          % ALIASING!!!
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% Doing it properly 
subsample_factor = 4;
 
tvec2 = decimate(tvec, subsample_factor);
y2 = decimate(y, subsample_factor);

subplot(212)          % No aliasing!!!
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% Neuralynx Data

cd /home/jalockli/Documents/BIOL680/data/R016-2012-10-08
fname = 'R016-2012-10-08-CSC03b.ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

%Don't expect to see spikes as the sampling rate is too slow to catch them.

% 16 bitts with a range of +_ 1500uV. Precision is:
1500 / ( 0.5 * (2^16) )

% So, the precision, or smallest voltage change that can be resolved would
% be about 0.05uV. ADBitVolts in hHeader is this value in Volts.

%From the header:
ADBitVolts = 4.57778e-008 

%Samples is in AD values, to convert to Volts:
% Assuming AD = 0 = Volts,
% 

Samples_V = Samples * ADBitVolts; % convert to Volts

Timestamps_s = Timestamps * 10^(-6); % convert to Seconds

%Total number of samples:
length(Samples(:))

%From the header: '-SamplingFrequency 2000' (us)
%Total duration assuming perfect sampling schedule:

sampling_Frequency = 2000 
length(Samples(:)) * sampling_Frequency

%Actual Time elapsed:
Timestamps_s(end) - Timestamps_s(1)

% Restricting data to rat running time
run(FindFile('*keys.m'))























