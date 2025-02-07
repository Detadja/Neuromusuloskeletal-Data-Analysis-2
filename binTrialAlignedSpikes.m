function [trial_spike_rate, time_bins] = binTrialAlignedSpikes(...
    aligned_spike_times, aligned_spike_labels, time_before, time_after, bin_width)

%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
% The function description below described the high-level goals of the function and formats of the inputs and outputs. Read this carefully.

%[trial_spike_rate, time_bins] = binTrialAlignedSpikes(aligned_spike_times, aligned_spike_labels, time_before, time_after, bin_width)
%
%function to bin trial-aligned spike times into trial-aligned spike rate
%estimates.
%A. Orsborn (last updated 12/19/20)
%
%inputs: aligned_spike_times - cell {#channels x 1} with vector (#spikes x 1)
%                               with trial-aligned time-stamps for each
%                               detected spike for a given unit
%         aligned_spike_labels - cell {#channels x 1} with vector (#spikes x 1)
%                                with label of the trial# each spike
%                                belongs to
%        time_before - time prior to alignment event to load for a
%                      given trial (in seconds)
%        time_after  - time after alignment event to load for a
%                      given trial (in seconds)
%        bin_width   - width of time-bin to use for estimating spike rate
%                      (in seconds)
%
%outputs: trial_spike_rate    - matrix [#trials x time x #channels] of estimated spike-rates for each trial/unit
%         time_bins           - vector (#time points x 1) of the center of
%                               the time-bin over which firing rate is estimated. (in seconds)
%

%get size of data
num_channels = length(aligned_spike_times);
num_trials = max( cat(1, aligned_spike_labels{:}));

%define the bin edges for your data
time_bins = -time_before : bin_width : (time_after);
num_time_points = length(time_bins);

%initialize a matrix for the binned spikes [# trials x time x #channels]
trial_spike_count = zeros(num_trials, num_time_points, num_channels); %FILLIN


%loop through channels
for iC = 1:num_channels

    %loop through trials
    for iT = 1:num_trials

        %find all the spikes for unit iC that belong to trial iT
        %hint: recall we saved a trial-label variable.
        spike_idx = (aligned_spike_labels{iC} == iT); %FILLIN - logical variable for does/doesn't belong to this trial
        data = aligned_spike_times{iC}(spike_idx);

        %bin the spike times in 'data' into counts using histc
        data_binned = histc(data, time_bins);

        %make sure data_binned is a row vector.
        %(If channel has 1 spike, histc returns a column vector)
        data_binned = data_binned(:);

        %deal with trials w/o any spikes (histc will return empty vector)
        if isempty(data_binned)
            data_binned = zeros(length(binvector),1);
        end

        %store binned spikes for current unit and trial into
        %trial_spike_count
        trial_spike_count(iT,:, iC) = data_binned; %fill in
    end

end

%remove the last bin (histc will return zero count for last bin always--meaningless)
trial_spike_count = trial_spike_count(:,1:end-1,:);
time_bins = time_bins(1:end-1) + diff(time_bins)/2; %change time_bin to be the center of the bin, not the edges.

%convert trial_spike_count to a rate (hertz)
trial_spike_rate = trial_spike_count / bin_width; %FILLIN

end
