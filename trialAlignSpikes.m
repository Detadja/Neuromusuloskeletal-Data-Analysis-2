function [aligned_spike_times, aligned_spike_labels] = ...
    trialAlignSpikes(spike_times, align_time, time_before, time_after)

%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
% The function description below described the high-level goals of the function and formats of the inputs and outputs. Read this carefully.

% [aligned_spike_times, aligned_spike_labels] = ...
%    trialAlignSpikes(spike_times, align_time, time_before, time_after)
%
%function to align spikes to trial-event-times
%
%A. Orsborn (last updated 12/19/20)
%
%inputs: spike_times - cell {#channels x 1} containing vector (#spikes x 1) of time-stamps
%                      for detected spikes for a given unit
%        align_time  - vector (#trials x 1) of time-stamps for trial events to align spikes to.
%        time_before - time prior to alignment event to load for a
%                      given trial (in seconds)
%        time_after  - time after alignment event to load for a
%                      given trial (in seconds)
%
%outputs: aligned_spike_times - cell {#channels x 1} with vector (#spikes x 1)
%                               with trial-aligned time-stamps for each
%                               detected spike for a given channel
%         aligned_spike_labels - cell {#channel x 1} with vector (#spikes x 1)
%                                with label of the trial# each spike
%                                belongs to
%


num_channels = length(spike_times);


%%%%%% trial-alignment for spikes

%find the alignment-event times
num_trials = length(align_time);

%loop through each spike variable
aligned_spike_times = cell(num_channels,1); %initialize variable for aligned spike-times
aligned_spike_labels = cell(num_channels,1); %initialize variable for trial-id for each spike (which trial each spike belongs to)
for iC=1:num_channels


    %loop through trials
    aligned_spike_times{iC} = []; %initialize empty vectors for spike-times and trial-ids
    aligned_spike_labels{iC} = [];
    for iT=1:num_trials

        %make a logical vector (spike_idx) that indicates all spike-times that occur
        %in the trial time-window: [time_before time_after] + align_time(iT)
        trial_start_time = align_time(iT) - time_before; %FILLIN
        trial_end_time = align_time(iT) + time_after; %FILLIN
        spike_idx = (spike_times{iC} >= trial_start_time & spike_times{iC} <= trial_end_time); %FILLIN

        %subtract the align_time(iT) from all spike times to set t=0 at
        %align_time for each trial
        tmp_aligned_spike_times = spike_times{iC}(spike_idx) - align_time(iT);

        %concatonate trial-aligned spikes
        aligned_spike_times{iC} = cat(1, aligned_spike_times{iC}, tmp_aligned_spike_times);

        %make a vector called tmp_id that is [#spikes in this trial x 1] = iT to keep track of
        %which trial each spike belongs to.
        tmp_id = repmat(iT, length(spike_idx(spike_idx == 1)), 1); %FILLIN
        aligned_spike_labels{iC} = cat(1, aligned_spike_labels{iC}, tmp_id); %make a master-list of labels by concatonating across trials

    end %loop through trials

end %loop through units
