% first, cd to where the data you just grabbed is located
cd /home/jalockli/BIOL680/data/R042-2013-08-018/

%% load the data (note, may need to unzip position data first)
fc = FindFiles('*.t');
S = LoadSpikes(fc);

[csc,csc_info] = LoadCSC('R042-2013-08-18-CSC03a.ncs');

[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );
 
%% plot


