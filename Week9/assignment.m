% Assignment Sandbox

% loading
cd /home/jalockli/Documents/BIOL680/data/R042-2013-08-18/
 
fc = FindFiles('*.t');
S = LoadSpikes(fc);


%% for each of the two neurons, restrict the data to [3200 5650] 
%(the time interval when the rat was running on the track)

cell1_id = 5; cell2_id = 42;
t = [3200 5650]; %seconds
%t = [5801 5801.7];

s1 = Data(Restrict(S{cell1_id},t(1),t(2))); 
s2 = Data(Restrict(S{cell2_id},t(1),t(2)));


%% compute the spike density function for each, 
%       making sure that your tvec runs from 3200 to 5650
%       ensure 50ms SD for the Gaussian convolution kernel

binsize = 0.001; % 1ms bins
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;

% s1
figure;
spk_count = histc(s1,tbin_edges);
spk_count = spk_count(1:end-1);

%Plot spikes
line([s1 s1],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
hold on;

% Plot Gaussian
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.05./binsize; % 50ms SD
gk = gausskernel(gauss_window, gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count, gk,'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');


% s2
figure;
spk_count = histc(s2,tbin_edges);
spk_count = spk_count(1:end-1);

%Plot spikes
line([s2 s2],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
hold on;

% Plot Gaussian
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.05./binsize; % 50ms SD 
gk = gausskernel(gauss_window, gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count, gk,'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');

 
%% to use these SDFs to generate Poisson spike trains, convert the firing rates given by the SDF to a probability of emitting a spike in a given bin. (As you did above for a 0.47 Hz constant firing rate.)
 
%% generate Poisson spike trains, making sure to use the same tvec
 
%% convert Poisson spike trains to ts objects and compute the ccf