function csc = LoadCSC(fname)
% function csc = LoadCSC(fname)
%
% improved version of LoadCSC() that handles missing data 
% (as flagged by NumberOfValidSamples) correctly,
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


TimeConvFactor = 10^-6; % from nlx units to seconds
VoltageConvFactor = 10^6; % from volts to microvolts
% extract_varargin;

fname
pwd
% load data
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

% extract information from header
csc_info = readCSCHeader(Header);




% Use NumberOfValidSamples to remove invalid samples
% First, a complete data block is:
full = max(NumberOfValidSamples) % is 512 in this example

% Create a vector of the number of zeros to trim from each block
bad_zeros = full - NumberOfValidSamples; %# zeros to erase in each column

% Indexes of zeros
index = zeros(size(Samples));
index( find(Samples == 0) ) = 1;   %An array locating all the zeros
totals = sum(index);        %Column totals of zeros (good and bad)
good_zeros = totals - bad_zeros; %# zeros to preserve in each column

%Go through only the columns with bad zeros, and replace them with NaNs
for column = find(bad_zeros)   
   [x y] = find(Samples(:,column) == 0, bad_zeros(column) ) ;
   Samples(x, column) = NaN ;
end

%Sanity check
index = zeros(size(Samples));
index( find(Samples == 0) ) = 1;
totals = sum(index);
if not( sum(totals == good_zeros) == length(totals) )
    printf('ERROR: There are problems with missing data')
end


% unwrap samples
csc_data = reshape(Samples,[size(Samples,1)*size(Samples,2) 1]);

% apply conversion
csc_data = csc_data.*VoltageConvFactor.*csc_info.ADBitVolts;

% create matching timestamps; idea is to make a matrix to match Samples
% original and then add dt's within each 512 sample packet

%Check for bad timestamp values, use interpolation to replace them:
while min(diff(Timestamps))< 0 %keep repeating until no more -ve diffs
    bad_timestamps = find(diff(Timestamps)< 0)
    for x = bad_timestamps
        printf('Interpolating to replace a bad timestamp');
        Timestamps(x) = mean(Timestamps(x-1), Timestamps(x+1));
    end
end



%pad the columns
csc_timestamps = repmat(Timestamps,[size(Samples,1) 1]).*TimeConvFactor;

%Create a matrix accounting for estimated time since timestamp
dtvec = (0:size(Samples,1)-1)*(1/csc_info.SamplingFrequency);
dtmat = repmat(dtvec',[1 size(Samples,2)]);

% add it to the timestamp matrix
csc_timestamps = csc_timestamps+dtmat;

%Now unwrap the timestamps just like the samples
csc_timestamps = reshape(csc_timestamps,[size(csc_timestamps,1)*size(csc_timestamps,2) 1]);



% Note that the header is included in the tsd object
%csc = mytsd(csc_timestamps, csc_data, csc_info);

% Using the standard tsd class: 
csc = tsd(csc_timestamps, csc_data);



function csc_info = readCSCHeader(Header)

csc_info = [];
for hline = 1:length(Header)
   
    line = strtrim(Header{hline});
    
    if isempty(line) | ~strcmp(line(1),'-') % not an informative line, skip
        continue;
    end
    
    a = regexp(line(2:end),'(?<key>\w+)\s+(?<val>\S+)','names');
    
    % deal with characters not allowed by MATLAB struct
    if strcmp(a.key,'DspFilterDelay_�s')
        a.key = 'DspFilterDelay_us';
    end
    
    csc_info = setfield(csc_info,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        csc_info = setfield(csc_info,a.key,str2double(a.val));
    end
    
end
