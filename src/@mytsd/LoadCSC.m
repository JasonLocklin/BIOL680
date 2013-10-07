% improved version of LoadCSC() that handles missing data 
% (as flagged by NumberOfValidSamples) correctly,

function csc = LoadCSC(fname)
% function csc = LoadCSC(fname)
%
% load a Neuralynx .ncs file correctly. Unlike the origional, it does not
% return timestamps with negative diffs, and no invalid samples.
%
% INPUTS:
% fname: [1 x n] char, i.e. a filename such as 'R042-2013-08-18-CSC01a.ncs'
%
% OUTPUTS:
%
% csc: [1 x 1] mytsd  %returns a mytsd object containing timestamps
% (seconds), corresponding samples (millivolts), and file header
% (csc_info). 



