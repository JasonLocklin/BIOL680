
function R = Range(tsa)

% tsd/Range
% 
%  R = Range(tsa)
%
%  returns range covered by tsa


% From sandbox:
% Use the first element of ExpKeys.TimeOnTrack and ExpKeys.TimeOffTrack 
% to find the indices of Timestamps corresponding to the Value session. 

[a RecordingStartIndex] = min(abs(Timestamps_s - ExpKeys.TimeOnTrack(1)))
[a RecordingEndIndex] =   min(abs(Timestamps_s - ExpKeys.TimeOffTrack(1)))


%Create a new set of variables on this restricted timeset
TimestampsValue = Timestamps_s(RecordingStartIndex:RecordingEndIndex);
SamplesValue    = Samples_V(1:512, RecordingStartIndex:RecordingEndIndex);
NumberOfValidSamplesValue = NumberOfValidSamples(RecordingStartIndex:RecordingEndIndex);

