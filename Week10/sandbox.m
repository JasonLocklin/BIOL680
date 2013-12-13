% Week 10 Tutorial

cd /home/jalockli/Documents/BIOL680/data/R042-2013-08-18/
 
sd.fc = FindFiles('*.t');
sd.fc = cat(1,sd.fc,FindFiles('*._t')); % also load poorly isolated cells (if you don't have these, get from database)
sd.S = LoadSpikes(sd.fc);

%% Neuralynx Video Data
% remember to unzip VT1.zip if not yet done
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, []);
Timestamps = Timestamps*10^-6;

% remove zero samples, which happen when no LED is detected
toRemove = (X == 0 & Y == 0);
X = X(~toRemove); Y = Y(~toRemove); Timestamps = Timestamps(~toRemove);

sd.x = tsd(Timestamps,X');
sd.y = tsd(Timestamps,Y');

% Using `Timestamps`,  determine the sampling rate of the video data.
min(diff(Timestamps))
mode(diff(Timestamps))
median(diff(Timestamps))

% Looks like the sampling rate isn't particularily precise or consistant.
% Since it agrees with the median value, a rounded modal rate is probably
% best (~33ms). I'm going to guess that the equipment is calibrated to a
% nice round frequency (rather than sampling interval), so it's probably 30Hz.

%% Visual Inspection
t = [3250 5650];
for iC = 1:length(sd.S)
    sd.S{iC} = Restrict(sd.S{iC},t(1),t(2));
end
sd.x = Restrict(sd.x,t(1),t(2));
sd.y = Restrict(sd.y,t(1),t(2));

plot(Data(sd.x),Data(sd.y),'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;

% Highlight a single cell (5)
% get x and y coordinate for times of spike
iC = 5;
spk_x = interp1(Range(sd.x),Data(sd.x),Data(sd.S{iC}),'linear');
spk_y = interp1(Range(sd.y),Data(sd.y),Data(sd.S{iC}),'linear');
 
h = plot(spk_x,spk_y,'.r');

%% Estimating tuning curves
SET_xmin = 10; SET_ymin = 10; SET_xmax = 640; SET_ymax = 480;
SET_nxBins = 63; SET_nyBins = 47;
 
spk_binned = ndHistc(cat(1,spk_x',spk_y'),[SET_nxBins; SET_nyBins],[SET_xmin; SET_ymin],[SET_xmax; SET_ymax]);
 
imagesc(spk_binned');
axis xy; colorbar;
