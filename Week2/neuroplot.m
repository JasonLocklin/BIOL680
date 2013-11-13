function neuroplot(spikes,csc,varargin)
% function neuroplot(spikes,csc,varargin)
%
% inputs:
%
% spikes: {nCells x 1} cell array of ts objects (spike train)
% csc: {nCSCs x 1} cell array of tsd objects (LFPs)
%
% varargins:
%
% cscColor: [nCSC x 3] RGB values to plot CSCs, default []
% spikeColor: [nCells x 3] RGB values to plot spikes, default []
%
% evt: [nEvents x 1] event times to plot, default []
% evtColor: [nEvents x 3] RGB values to plot events, default []
%
% interactiveMode: boolean to enable/disable arrow key navigation

cscColor   = [0.2 0.2 0.8]; %Set defaults for optional arguments
spikeColor = [0.4 0.4 0.7];
evt        = [];
evtColor   = [0.8 0.2 0.2];
interactiveMode = true;
%default x range is 2 seconds

extract_varargin; % override optional arguments

if   interactiveMode, f = figure('KeyPressFcn', @navigate);
else                  f = figure();
end

%Turn the data into two vectors (X and Y) for ploting.
spiked = [];
yvalues = [];
for yindex = 1:length(spikes)
    spiked = vertcat(spiked, Data(spikes{yindex}));
    yvalues = vertcat(yvalues, ones(length(Data(spikes{yindex})),1)*yindex);
end

% Scatterplot is very fast, but only has a small set of markers available.
% In another language, I would use |, but here, I'll use a triangle.
scatter(spiked, yvalues, 25, spikeColor, '>');
set(gca, 'XLim', [min(spiked) min(spiked)+2]); %Start with the first 2 sec
set(gca, 'ytick', []); %Hide Y axis

%Calculate scaling factor for LFPs (max range +-10)
cf = diff(minmax(Data(csc)))/20;
bias = get(gca, 'YLim');
bias = bias(2) + 10;

set(gca, 'YLim', get(gca, 'YLim')+ [0 20]); %Make space on the plot
hold on;
pl = plot(Range(csc), (Data(csc)/cf)+bias, 'Color', cscColor);

for i = evt
    line([2140 2140], get(gca, 'YLim'), 'Color', evtColor)
end

