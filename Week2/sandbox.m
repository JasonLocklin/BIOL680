set(0,'DefaultFigureWindowStyle','docked');

% first, cd to where the data you just grabbed is located
cd /home/jalockli/BIOL680/data/R042-2013-08-18/

%% load the data (note, may need to unzip position data first)
fc = FindFiles('*.t');
S = LoadSpikes(fc);

[csc,csc_info] = LoadCSC('R042-2013-08-18-CSC03a.ncs');

[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );

%%   ts and tsd objects: basic methods
% verify csc object has timestamps and data
if length(Data(csc)) == length(Data(csc));
    fprintf('csc is good\n')
else
    fprintf('csc is bad\n')
end

minmax(Range(csc))
csc_small = Restrict(csc, 5950, 6050);
length(Data(csc_small)) 

x = tsd((Timestamps*10^-6), X');
y = tsd((Timestamps*10^-6), Y');
x_small = Restrict(x, 5950, 6050);
length(Data(x_small))

%% Plotting with handles and figure callback function (navigate.m)
f = figure('KeyPressFcn', @navigate);
pl = plot(Range(csc_small), Data(csc_small));

set(gcf,'Color',[0 0 0]);
set(gca, 'XColor', [1,1,1]);
set(gca, 'YColor', [1,1,1]);
set(gca, 'XLim', [5989, 5990], 'FontSize', 24);

hold on; box off;

csc_mean = nanmean(Data(csc));
xr = get(gca, 'XLim');
mean_hdl = plot(xr, [csc_mean csc_mean]);
set(mean_hdl, 'Color', [1,0,0]);
set(mean_hdl, 'LineStyle', '--');
set(mean_hdl, 'LineWidth', 2);

%% Print it
set(gcf, 'InvertHardcopy', 'off')
print(gcf,'-dpng','-r300','R042-2013-08-18-LFPsnippet.png');

%% Anonymous Functions
sqr_fn = @(x) x.^2;
sqr_fn(2)

%% Neuroplot
% S is spikes

%if interactive:
f = figure('KeyPressFcn', @navigate);

spikes = [];
yvalues = [];

%Turn the data into two vectors (X and Y) for ploting.
for yindex = 1:length(S)
    spikes = vertcat(spikes, Data(S{yindex}));
    yvalues = vertcat(yvalues, ones(length(Data(S{yindex})),1)*yindex );
end

% Scatterplot is very fast, but only has a small set of markers available
% In another language, I would use |, but here, I'll just use a triangle.
scatter(spikes, yvalues, '>');
set(gca, 'XLim', [min(spikes) min(spikes)+2]);
%set(gca,'box','off','ycolor',[0.5 0.5 0.5]);
% set(gca, 'ytick', []);


%LFPs
hold on;

%Calculate scaling factor for max range +-10
cf = diff(minmax(Data(csc)))/20;
bias = get(gca, 'YLim');
bias = bias(2) + 10;

set(gca, 'YLim', get(gca, 'YLim')+ [0 20]);
pl = plot(Range(csc), (Data(csc)/cf)+bias );





