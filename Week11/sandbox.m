% Week 11 Tutorial

% In the shell:
% > cd /home/jalockli/Documents/BIOL680/data/
% > cp -r R016-2012-10-03/ R016-2012-10-03_tmp
% > cd R016-2012-10-03_tmp/
% > rename 's/\._t$/.t/' *

% cd to R016-2012-10-03 folder
cd /home/jalockli/Documents/BIOL680/data/R016-2012-10-03_tmp/
 
sd.fc = FindFiles('*.t');
sd.S = LoadSpikes(sd.fc);
 
csc = LoadCSC('R016-2012-10-03-CSC02d.ncs');
 
cscR = Range(csc); cscD = Data(csc);
 
Fs = 2000; dt = 1./Fs;
cscD = locdetrend(cscD,Fs,[1 0.5]); % remove slow drifts in signal (this can mess up the STA)


%% grab a LFP snippet for each spike and average (SLOW!!!)
w = [-1 1]; % time window to compute STA over
tvec = w(1):dt:w(2); % time axis for STA
 
iC = 3; % only do the third neuron for now
clear sta;
 
spk_t = Data(sd.S{iC});
 
h = waitbar(0,sprintf('Cell %d/%d...',iC,length(sd.S)));
 
for iSpk = length(spk_t):-1:1 % for each spike...
 
   sta_t = spk_t(iSpk)+w(1);
   sta_idx = nearest(cscR,sta_t); % find index of leading window edge
 
   toAdd = cscD(sta_idx:sta_idx+length(tvec)-1); % grab LFP snippet for this window
   % note this way can be dangerous if there are gaps in the data
 
   sta{iC}(iSpk,:) = toAdd'; % build up matrix of [spikes x samples] to average later
 
   waitbar(iSpk/length(spk_t));
end
 
close(h);

%% Now plot it
plot(tvec,nanmean(sta{3}),'k','LineWidth',2); 
set(gca,'FontSize',14,'XLim',[-0.5 0.5]); xlabel('time (s)'); grid on;

% There seems to be a certain asymmetry in the STA theta oscillation 
% (clearer on the left than on the right). How might we interpret this? 

% Signal to the left of zero in the STA plot indicate characteristics of
% the LFP prior to a typical spike, while signal to the left indicates
% properties of the LFP after a typical spike. It seems that, in this case,
% there is a more reliable theta oscilation prior to a typical spike. This
% may provide some evidence that theta oscilations are causal factors in
% spikes of this particular neuron, rather than the reverse.

%% Profiler
profile on

%% Profiler off
profile off

%% Faster STA computation that doesn't require a loop
bin_edges = cscR+(dt/2);
len = length(tvec);
 
clear sta;
iC = 1;
 
spk_ts = Restrict(sd.S{iC},cscR(1)-w(1),cscR(end)-w(2));
spk_t = Data(spk_ts)+w(1); % times corresponding to start of window
 
[~,spk_bins] = histc(spk_t,bin_edges); % index into data for start of window
 
spk_bins2 = spk_bins(:, ones(len,1));
toadd = repmat(0:len-1,[length(spk_bins) 1]);
 
spk_bins3 = spk_bins2+toadd;
 
sta = cscD(spk_bins3);






















