%% Week 8 Tutorial


%% Are the two signals in the lower left panel coherent? Ignore the amplitude changes for now. 
% Yes, the peaks of the upper wave tend to coincide with the point just
% after the troughs of the lower wave.

%% Diversion: the Wiener-Khinchin theorem
Fs = 500; dt = 1./Fs;
t = [0 2]; tvec = t(1):dt:t(2)-dt;
 
f1 = 8;
data1 = sin(2*pi*f1*tvec)+0.1*randn(size(tvec));
 
[acf,lags] = xcorr(data1,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
plot(lags,acf); grid on;

%% Cross-correlation
f2 = 8;
data2 = sin(2*pi*f2*tvec+pi/4)+0.1*randn(size(tvec)); % phase-shifted version of data1
 
[ccf,lags] = xcorr(data1,data2,100,'coeff'); % now a cross-correlation
lags = lags.*(1./Fs); % convert samples to time
plot(lags,ccf); grid on;
 %Note peak is offset from zero by phase difference between signals
 
%% Verify that changing the phase shift in data2 indeed changes the phase of the cross-correlogram. 
f2 = 8;
shift = pi
data3 = sin(2*pi*f2*tvec+shift)+0.1*randn(size(tvec)); % phase-shifted version of data1
 
[ccf,lags] = xcorr(data1,data3,100,'coeff'); % now a cross-correlation
lags = lags.*(1./Fs); % convert samples to time
plot(lags,ccf); grid on;
 %Note that the signals are completely out of phase, a trough is at zero.
 
%% Definition and example

figure;
subplot(221);
plot(tvec,data1,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1,hanning(250),125,length(data1),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1),Fs);
subplot(222)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');
 
[Pxy,F] = cpsd(data1,data2,hanning(250),125,length(data1),Fs);
subplot(223)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');
 
[acf,lags] = xcorr(data1,data2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(224)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');

%% What happens if you change the amplitude of one of the input signals? Make the amplitude of data1 twice as large. 
figure;
subplot(221);
plot(tvec,data1*2,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1*2,hanning(250),125,length(data1*2),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1*2),Fs);
subplot(222)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');
 
[Pxy,F] = cpsd(data1*2,data2,hanning(250),125,length(data1*2),Fs);
subplot(223)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');
 
[acf,lags] = xcorr(data1*2,data2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(224)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');

% Answer: The power level of the PSD is 10000x greater, the
% x-spectrum is 100x, and xcorr is identical.

%% For the two signals of unequal amplitude, compute the coherence by normalizing the cross-spectrum

figure;
subplot(241);
plot(tvec,data1,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1,hanning(250),125,length(data1),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1),Fs);
subplot(242)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');

[Pxy,F] = cpsd(data1,data2,hanning(250),125,length(data1),Fs);
subplot(243)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');


C = (abs(Pxy).^2)./(Pxx.*Pyy);
subplot(244)
plot(F,C); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('coherence');
 
% now at double amplitude
subplot(245);
plot(tvec,data1*2,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1*2,hanning(250),125,length(data1*2),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1*2),Fs);
subplot(246)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');

[Pxy,F] = cpsd(data1*2,data2,hanning(250),125,length(data1*2),Fs);
subplot(247)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');


C = (abs(Pxy).^2)./(Pxx.*Pyy);
subplot(248)
plot(F,C); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('coherence');


%% Find out how small you need to make data1 relative to the noise before a drop in coherence occurs. 
figure;
subplot(241);
plot(tvec,data1,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1,hanning(250),125,length(data1),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1),Fs);
subplot(242)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');

[Pxy,F] = cpsd(data1,data2,hanning(250),125,length(data1),Fs);
subplot(243)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');


C = (abs(Pxy).^2)./(Pxx.*Pyy);
subplot(244)
plot(F,C); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('coherence');
 
% now shrink data1
subplot(245);
factor = 10.^-245 % I don't hink this is right. I think we are just moving outside the range of the FP math.
plot(tvec,data1.*factor,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1.*factor,hanning(250),125,length(data1.*factor),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1.*factor),Fs);
subplot(246)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');

[Pxy,F] = cpsd(data1.*factor,data2,hanning(250),125,length(data1.*factor),Fs);
subplot(247)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');

%Calculate coherence
C = (abs(Pxy).^2)./(Pxx.*Pyy);
subplot(248)
plot(F,C); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('coherence');


%% just verify some cases where we break the phase relationship
f = 2; % freq modulation (Hz) 
f2 = 8;
m = 2; % freq modulation strength
wsz = 200; % window size 
 
subplot(421)
s2 = data2;
plot(tvec,s2,tvec,data1); title('signal 1 - constant phase');
 
subplot(422)
s3 = sin(2*pi*f2*tvec + m.*sin(2*pi*f*tvec - pi/2)) + 0.1*randn(size(tvec));
plot(tvec,s3,tvec,data1); title('signal 2 - varying phase');
 
subplot(423)
[Ps2,F] = pwelch(s2,hanning(wsz),wsz/2,length(data2),Fs);
plot(F,abs(Ps2)); title('PSD');
 
subplot(424)
[Ps3,F] = pwelch(s3,hanning(wsz),wsz/2,length(data2),Fs);
plot(F,abs(Ps3)); title('PSD');
 
subplot(425)
[C,F] = mscohere(data1,s2,hanning(wsz),wsz/2,length(data1),Fs); % shortcut to obtain coherence
plot(F,C); title('coherence'); xlabel('Frequency (Hz)');
 
subplot(426)
[C,F] = mscohere(data1,s3,hanning(wsz),wsz/2,length(data1),Fs);
plot(F,C); title('coherence'); xlabel('Frequency (Hz)');
 
[acf,lags] = xcorr(data1,s2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(427)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');
 
[acf,lags] = xcorr(data1,s3,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(428)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');

%% What happens when you change the window size to 500ms? 
% cross-correlations are sharper/higher, but the coherence signal is lost.
% in noise.


%% Overall comparison of vStr-vStr and vStr-HC coherence
%cd R016-2012-10-03/ %?

run(FindFile('*keys.m'));
vStr1_csc = LoadCSC(ExpKeys.goodGamma{1});
vStr2_csc = LoadCSC(ExpKeys.goodGamma{2});
HC_csc = LoadCSC(ExpKeys.goodTheta{1});
 
% restrict to shorter segment for speed
vStr1_csc = Restrict(vStr1_csc,ExpKeys.TimeOnTrack(2),ExpKeys.TimeOffTrack(2));
vStr2_csc = Restrict(vStr2_csc,ExpKeys.TimeOnTrack(2),ExpKeys.TimeOffTrack(2));
HC_csc = Restrict(HC_csc,ExpKeys.TimeOnTrack(2),ExpKeys.TimeOffTrack(2));
 
vStr1_csc = Data(vStr1_csc); vStr2_csc = Data(vStr2_csc); HC_csc = Data(HC_csc);
%%

Fs = 2000; wsize = 2048;
 
% compute PSDs and coherences for each signal and each pair respectively
[P1,F] = pwelch(vStr1_csc,hanning(wsize),wsize/2,2*wsize,Fs);
[P2,F] = pwelch(vStr2_csc,hanning(wsize),wsize/2,2*wsize,Fs);
[P3,F] = pwelch(HC_csc,hanning(wsize),wsize/2,2*wsize,Fs);
[C1,F] = mscohere(vStr1_csc,vStr2_csc,hanning(wsize),wsize/2,2*wsize,Fs);
[C2,F] = mscohere(vStr1_csc,HC_csc,hanning(wsize),wsize/2,2*wsize,Fs);
 
% plot
subplot(121)
h(1) = plot(F,10*log10(P1),'k','LineWidth',2); hold on;
h(2) = plot(F,10*log10(P2),'g','LineWidth',2); hold on;
h(3) = plot(F,10*log10(P3),'m','LineWidth',2); hold on;
set(gca,'XLim',[0 150],'XTick',0:25:150,'FontSize',12); grid on;
legend(h,{'vStr1','vStr2','HC'},'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 
 
subplot(122)
h(1) = plot(F,C1,'LineWidth',2); hold on;
h(2) = plot(F,C2,'r','LineWidth',2);
set(gca,'XLim',[0 150],'XTick',0:25:150,'FontSize',12); grid on;
legend(h,{'vStr1-vStr2','vStr1-HC'},'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Coherence');


%% Comparison of vStr-HC coherence between experimental conditions
cd /home/jalockli/Documents/BIOL680/data/R016-2012-10-03/
fc = {'R016-2012-10-03-CSC04d.ncs','R016-2012-10-03-CSC03d.ncs','R016-2012-10-03-CSC02b.ncs'};

% Note: in order for ft_read_neuralynx to work, I hadd to add:
%       key = regexprep(key,'[^a-zA-Z]','');
% to the 77th line of neuarlynx_getheader to remove unicode characters
% from the headers. My version of matlab does not allow unicode characters
% for variable identifiers (including as struct keys in this case).
% Even though it allows unicode in the term, editor, and string contents.

data = ft_read_neuralynx_interp(fc);
data.label = {'vStr1','vStr2','HC1'}; % reassign labels to be more informative, the original filenames can be retrieved from the header


cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5; cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk', 'both'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);

%% compute trial-averaged cross-spectrum
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2' 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = -2:0.05:0; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_pre = ft_freqanalysis(cfg, data_trl);

%% Compute coherence from cross-spectrum
cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_pre);

figure;
cols = 'rgb';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
 
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
legend(h,lbl);

%% Do the same for the 0-2 second post-nosepoke period

cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2' 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = 0:0.05:2; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_post = ft_freqanalysis(cfg, data_trl);

cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_post);

figure;
cols = 'rgb';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
 
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
legend(h,lbl);


%% Comparing the coherence spectra for 1 and 5 pellet trials

% 1 Food Pellet
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5; cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk', 'both'
cfg.trialdef.cue = {'c1'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);
% X-sepectrum
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2' 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = 0:0.05:2; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_pre = ft_freqanalysis(cfg, data_trl);
% Coherence
cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_pre);

figure;
cols = 'rgb';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
 
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
legend(h,lbl);

% 5 Food Pellets
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5; cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk', 'both'
cfg.trialdef.cue = {'c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);
% X-sepectrum
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2' 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = 0:0.05:2; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_pre = ft_freqanalysis(cfg, data_trl);
% Coherence
cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_pre);

figure;
cols = 'rgb';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
 
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
legend(h,lbl);

% There is a definite 30Hz spike in the post-1-pellet trials. For the 5-pelet
% trials, there is no real solid spike like the 1-pellet trials -perhapse 
% still a bit of a 30Hz spike along with a 41Hz spike, but I would be 
% cautious to say anything about it without more examination

%% Time-frequency coherence 
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2' 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = 0:0.05:2; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_post = ft_freqanalysis(cfg, data_trl);

% Compute coherence from cross-spectrum
cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_post);
iC = 1; % which signal pair to plot
lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.cohspctrm(iC,:,:))); axis xy; colorbar
xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);


%% Fixed time-window:

cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = repmat(1,1,100);  % Time dependent: 1sec time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2' 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = 0:0.05:2; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_post = ft_freqanalysis(cfg, data_trl);

% Compute coherence from cross-spectrum
cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_post);

iC = 1; % which signal pair to plot
lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.cohspctrm(iC,:,:))); axis xy; colorbar
xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);

% Certainly doesn't have the temporal change from low to high frequencies, 
% But there is still lots of seemingly spurious noise.

%% Beyond Coherence
cfg            = [];
cfg.method     = 'ppc';
fd             = ft_connectivityanalysis(cfg,TFR_post);
iC = 1; % which signal pair to plot

lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.ppcspctrm(iC,:,:))); axis xy; colorbar
xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);

% Much less spurious looking noise, and perhapse a 52Hz coherence signal
% lasting about 0.3 seconds, centered at 0.25s post stimuli.





