% Week 9 Sandbox file

%% Binned firing rates
% loading
cd /home/jalockli/Documents/BIOL680/data/R042-2013-08-18/
 
fc = FindFiles('*.t');
S = LoadSpikes(fc);
 
% plot
iC = 47;
t = [5801 5801.7];
 
spk_t = Data(Restrict(S{iC},t(1),t(2))); % get spike times
 
line([spk_t spk_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command

%% Bin the spikes
binsize = 0.1;                   % Size of bins
tbin_edges = t(1):binsize:t(2); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)
 
spk_count = histc(spk_t,tbin_edges); % get spike counts for each bin
spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.

hold on;
h = bar(tbin_centers,spk_count);
set(h,'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]); % reformat bar appearance
 
yt = get(gca,'YTick'); ytl = get(gca,'YTickLabel');
set(gca,'YLim',[-1.5 10],'YTick',yt(2:end),'YTickLabel',ytl(2:end,:)); xlabel('time (s)');


% The number of spikes per bin seems to align with the height of the
% histogram bars by visual inspection.

%% Reduce the bin size to 10 ms and replot. 

binsize = 0.01;                   % Size of bins
tbin_edges = t(1):binsize:t(2); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)
 
spk_count = histc(spk_t,tbin_edges); % get spike counts for each bin
spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.

hold on;
h = bar(tbin_centers,spk_count);
set(h,'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]); % reformat bar appearance
 
yt = get(gca,'YTick'); ytl = get(gca,'YTickLabel');
set(gca,'YLim',[-1.5 10],'YTick',yt(2:end),'YTickLabel',ytl(2:end,:)); xlabel('time (s)');

% The histogram now gives a decent representation of where the spikes are.


%% 1 ms and replot. 

binsize = 0.001;                   % Size of bins
tbin_edges = t(1):binsize:t(2); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)
 
spk_count = histc(spk_t,tbin_edges); % get spike counts for each bin
spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.

hold on;
h = bar(tbin_centers,spk_count);
set(h,'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]); % reformat bar appearance
 
yt = get(gca,'YTick'); ytl = get(gca,'YTickLabel');
set(gca,'YLim',[-1.5 10],'YTick',yt(2:end),'YTickLabel',ytl(2:end,:)); xlabel('time (s)');

%This basically gives us the spike train back.


%% Spike density functions (Convolution)
binsize = 0.001; % select a small bin size for good time resolution
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
spk_count = histc(spk_t,tbin_edges);
spk_count = spk_count(1:end-1);
 
rw_sdf = conv2(spk_count,rectwin(50),'same'); % convolve with rectangular window
%The same option specified in the call to conv2() ensures that the output is the same size as the input.
plot(tbin_centers,rw_sdf,'b');


%% Gaussian Window:

gau_sdf = conv2(spk_count,gausswin(50),'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');
% Note: convolved spike trains are referred to as spike density functions
% (SDF)

%% Gausskernel to scale to spikes per second, and specifying convolution window and Gaussian SD.
binsize = 0.001; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.02./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');

%% ISI histograms and the coefficient of variation

iC = 47;
spk_t = Data(S{iC}); % spike times
isi = diff(spk_t); % interspike intervals
 
dt = 0.001; % in s, because spike times are in s
isi_edges = 0:dt:0.25; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
 
isih = histc(isi,isi_edges);
 
bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;


%% Plot each ISI against it's successor (ISI return plot)
isi_shifted = zeros(size(isi));
isi_shifted(1:end-1) = isi(2:end);
isi_shifted(end) = isi(1);

figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
xlim(axes1,[0 0.25]);
ylim(axes1,[0 0.25]);
box(axes1,'on');
hold(axes1,'all');
plot(isi, isi_shifted, '.')
% Note: a log-log plot might be more useful to visualise where the bulk of
% the data lies.

%% Random synthetic spike train with Poisson point process
dt = 0.001;
t = [0 10]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);
 
pspike = 0.5; % probability of generating a spike in bin
rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
 
line([spk_poiss_t spk_poiss_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
axis([0 0.1 -1.5 5]); set(gca,'YTick','');

%% Plot the ISI histogram for synthetic spike train
iC = 47;
% spike times from spk_poiss_t
isi = diff(spk_poiss_t); % interspike intervals
 
dt = 0.001; % in s, because spike times are in s
isi_edges = 0:dt:0.25; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
 
isih = histc(isi,isi_edges);
 
bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.02]); xlabel('ISI (s)'); ylabel('count'); grid on;

% What should the probability of a spike be, for each 1ms bin, for our 
% synthetic spike train to have a mean firing rate of 0.47 spikes/s?

% A mean firing rate of 0.47 spikes/s is 0.00047 spikes per ms on average.
% The probability of a spike in a given 1ms bin is then 0.00047, or 0.047%


%% Recompute a spike train generated with this probability, instead of 0.5
%note: need to increase the length of the spike train to 4500s to get a number of spikes similar to that in the real spike train. 

dt = 0.001;
t = [0 4500]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

pspike = 0.00047; % probability of generating a spike in bin
rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time

%Real spike train:
1./mean(diff(spk_t))
%Synthetic:
1./mean(diff(spk_poiss_t))

%They are close. If time interval is bumped up higher, it get's closer.


%% Plot ISI return plot
isi = diff(spk_poiss_t)
isi_shifted = zeros(size(isi));
isi_shifted(1:end-1) = isi(2:end);
isi_shifted(end) = isi(1);

figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
xlim(axes1,[0 0.25]);
ylim(axes1,[0 0.25]);
box(axes1,'on');
hold(axes1,'all');
plot(isi, isi_shifted, '.')

% Note: no structure.

%% Verify that the synthetic, but not the real spike train has a refractory period
min(diff(spk_t))
min(diff(spk_poiss_t))

% The real spike train has a refractory period in the ballpark of 3ms,
% while the synthetic spike train has spikes less than 1ms apart.



%% Create an inhomogeneous Poisson spike train 

dt = 0.001;
t = [0 4500]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

% Create a vector of length tvec, oscilating between 0 and 0.5
pspike = ( sin( (1:length(tvec))/ (2*pi) ) + 1 ) .* 0.25;

rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time

%Plot a selection of the spikes to verify oscilating density
figure;
line([spk_poiss_t(1:20) spk_poiss_t(1:20)],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command


%% Spike autocoorrelation function

% Generate a poisson spike train with:
%        mean firing rate = 0.47 HZ
%        time bin         = 10 ms, from -1 to 1s


dt = 0.01;  %10ms time bin
t = [-1 1]; % 2 second interval (-1 to 1s)
tvec = t(1):dt:t(2);

pspike = 0.47; % probability of generating a spike in bin
%rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time  

% Check
figure;
line([spk_poiss_t spk_poiss_t],[-1 1],'Color',[0 0 0]); % note, plots all spikes in one command
length(spk_poiss_t) ./2 % spkes / second
 
% Convert to ts object
spk_poiss_ts = ts(spk_poiss_t)

% apply acf.m
[ac,xbin] = acf(spk_poiss_ts, 0.01, 1)
plot(xbin,ac)

% I don't get a flat autocorrelation like what is on the wiki. Not sure
% why.

%% Actual data
[acorr,xbin] = acf(S{47},0.01,1);
plot(xbin, acorr)
% with actual data, there appears to be some periodicity to the data, with
% significan lobes at around +-0.15s.


%% Spike cross-correlation function (ccf)
% What do you expect the cross-correlation function of two Poisson spike trains to look like?

% Answer: Two spike trains that "look alike" are two trains that have a
% similar pattern, but may be shifted somewhat relative to each other. The
% resulting x-correlation would therefore have a strong central peak
% somewhere in the vascinity of, but not necissarily at, zero.

% We can test this with our inhomogeneous Poisson spike train

dt = 0.01;
t = [-1 1]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);

% Create a vector of length tvec, oscilating between 0 and 0.5
pspike = ( sin( (1:length(tvec))/ (2*pi) ) + 1 ) .* 0.25;

rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t1 = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
spk_poiss_ts1 = ts(spk_poiss_t1)

spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t2 = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
spk_poiss_ts2 = ts(spk_poiss_t2)

%Verify the two visually simmilar poisson spike trains:
figure;
line([spk_poiss_t1 spk_poiss_t1],[-1 -0.75],'Color',[1 0 0]); % note, plots all spikes in one command
line([spk_poiss_t2 spk_poiss_t2],[-0.75 -0.5],'Color',[0 0 1]); % note, plots all spikes in one command


figure;
[cc, xbin] = ccf(spk_poiss_ts1, spk_poiss_ts2, 0.01, 1); 
plot(xbin, cc);

% Because there is no shift in the pattern relative to one another, there
% is a strong lobe centered around zero. Because of the oscilating pattern
% common to the two of them, there are lobes approximately +-310 and
% +-620ms.

%% Xcorr of two hippocampal place cells
cell1_id = 5; cell2_id = 42;
 
s1 = Restrict(S{cell1_id},3200,5650); % restrict to on-track times only
s2 = Restrict(S{cell2_id},3200,5650);
 
[xcorr,xbin] = ccf(s1,s2,0.01,1);
 
plot(xbin,xcorr);
set(gca,'FontSize',20); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('%d-%d',cell1_id,cell2_id));

