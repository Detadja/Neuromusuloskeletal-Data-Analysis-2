

function [spike_times, spike_waveforms] = detectSpikes(data, threshold, Fs)

  %%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
  % The function description below described the high-level goals of the function and formats of the inputs and outputs. Read this carefully.

% [spike_times, spike_waveforms] = detectSpikes(data, threshold, FS)
% Function to detect spikes in an electrophysiology trace by finding times when data amplitude (absolute value) exceeds a threshold.
% inputs: data - vector time-series of voltage input ([1 x time] or [time x 1])
%         threshold - scalar threshold value for spike detection (units same as data)
%         Fs        - sampling rate (samples/second) of 'data' vector
% outputs: spike_times - vector ([# spikes x 1]) of times (in seconds) of detected spikes
%          spike_waveforms - matrix ([# spikes x 30]) of waveform of each detected spike.
%


%define constants for the function
num_waveform_points_total = 30; %# of samples to load to capture the waveforms
num_waveform_points_before = 10; %# of samples before detection point to load

% we want to detect the time when the signal amplitude goes from below to above our threshold.
%1. create a logical vector that is zero when data is smaller than the threshold and 1 when it is bigger than the threshold
%we will want to account for the fact that we may want to use a negative threshold,
%in which case we want to find when the signal goes BELOW the threshold.
if threshold > 0
    amp_above_threshold = data > threshold; %FILLIN
else
    amp_above_threshold = data < threshold; %FILLIN
end

%2. create a vector that = 1 when amp_above_threshold goes from 0 -> 1 (hint: can be done in one line)
Low_hi = diff(amp_above_threshold); %ILLIN

%3. get the time index when low_hi = 1 (i.e. when our signal goes from below to above the threshold)
% this corresponds to the index of detected spikes
spike_idx = find(Low_hi); %FILLIN index in data
spike_times = spike_idx / Fs; %FILLIN convert to time in seconds

%now that we have detected spikes, we want to get the data around the time of a spike--the waveforms
num_spikes = length(spike_times); %# of spikes detected
spike_waveforms = nan(num_spikes, num_waveform_points_total);  %initial matrix for the waveforms

%loop through spikes
for iSpike=1:num_spikes
    %calculate the indices of data to load
    idx = spike_idx(iSpike) + [0:num_waveform_points_total-1] - num_waveform_points_before ;

    %check that idx is within bounds of data length (at edges may run out of data)
    if ~any( idx<0) & ~any(idx > length(data) )
        %save the calculated portion of 'data' into spike_waveforms
        spike_waveforms(iSpike,:) = data(idx); % FILL IN
    end
end %end loop through num_spikes

end %end function
