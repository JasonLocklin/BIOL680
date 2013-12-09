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


%% Plot each ISI against it's successor
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

%% Recompute the ISI histogram for a spike train generated with this probability, instead of 0.5
%note: need to increase the length of the spike train to 4500s to get a number of spikes similar to that in the real spike train. 




