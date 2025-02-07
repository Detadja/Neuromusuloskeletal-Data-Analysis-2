%%%%% This script provides a guided outline for analyzing your experimental data collected for Experiment 2 (quantifying sensory response to tactile stimulation)
%%%%% Written by A.L. Orsborn, v201219
%%%%%
%%%%%
%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.


% we will first load a data file and test our pre-processing and calculations on one file. Once that is complete, we can extend our analysis to all data to examine trends.

%% testing pre-processing and calculations

% Define some basic things to make it easy to find your data files.
% We will want to take advantage of systematic naming structure in our data files.
% Your files should have names like [prefix][date][id #].
% Note that the SpikeRecorder program automatically saves files with date and time in the name.
% We recommend re-naming your files to convert time into a simpler id# e.g. 1, 2, 3...

dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 2\Lab 2 (Data)\'; %FILLIN in with the path to where your data is stored

file_prefix = '2.'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

% file_date   = ''; %FILLIN with the date string used in your file
file_idnum  = 2; %FILLIN with the numeric value of the file id

full_file_name = [dataDir file_prefix num2str(file_idnum) file_type];


% load your data file and the sampling rate into matlab into the variables 'data' and FS, respectively.
[data, FS] = audioread(full_file_name);

% check the basic properties of your loaded variables
whos data FS

%note that if you are recording 2 channels instead of 1, your data will have 2 columns (#time points x #channels). Let's create a variable to keep track of how many channels we have so our code is flexible:
num_channels = 2; %FILLIN

% plot 1 second of your data (channels 1 and 2) to assure the loaded file looks correct
figure
time = [0:1/FS:1];
plot( time, data(1:length(time)) ) %fill in the portion in the brackets to only plot 1 second (hint: you need to take the sampling rate into account)
xlabel('Time, in samples') %fill in the units of the time axis on this plot
ylabel('Voltage (mV)')


% Now let's find spikes and their waveforms in the data.
% use the code you developed in lab 1 for this task (detectSpikes.m)
%note that you may need to modify this function or how you call it if you have two channels of data in your file.
%create a variable spike_times - cell {# channels x 1} with the detected spike times for each channel.
%note: if you recorded multiple channels, you might need a for loop to create the spike_times cell.
%FILLIN with relevant code
spike_times_cell = cell(num_channels, 1);
spike_waveforms_cell = cell(num_channels, 1);
channels = [1 2];
time_waveforms = [-10:19]./FS;

threshold = [-0.08 -0.2];
for iF=1:num_channels
    [spike_times_cell{iF}, spike_waveforms_cell{iF}] = detectSpikes(data(:, iF), threshold(iF), FS);

    figure(iF);
    plot(time_waveforms, spike_waveforms_cell{iF}) %FILLIN with relevant variable name
    xlabel('Time (s)') %FILL IN label
    ylabel('Waveforms (mV)') %FILL IN label, be sure to include units
    set(gca, 'ylim') %keep axis the same across plots for easier comparison
    title(['Channel ', num2str(channels(iF))])
    grid on
    % figure
    % time_data = [0:length(data(:,1))-1]./FS; %this creates a time vector of the data in seconds
    % plot(time_data, data(:,iF)) %plot the raw data
    % xlabel('Time (s)')
    % ylabel('Voltage (mV)')
    % hold on
    % plot([spike_times_cell{iF} spike_times_cell{iF}]', repmat([.1 .15], length(spike_times_cell{iF}),1)', 'r')
    % grid on
end

%inspect your data to convince yourself the function is working as you expect
%Some of the plots we generated in lab 1 may be helpful here (e.g. raw
%traces + detected spike times)



%%%%%
%now we need to learn how to load the 'events' files and analyze our spike data with respect to these events.

%1. load an 'events' file.
%We have provided a shell function for this operation called 'readEventsFile.m' Open it and fill it in.
events_file_type = '.txt';
events_file_name = [dataDir file_prefix num2str(file_idnum) events_file_type];

[EVENTS, EVENT_TIMES] = readEventsFile(events_file_name);


%2. We now need to write code to "align" our spike times with stimulation trials and make raster plots.
%we've provided a shell function for this operation called
%'trialAlignSpikes.m'. Open this function, fill it in, and then use it here.
align_time = EVENT_TIMES;
time_before = 1;
time_after  = 1;
[aligned_spike_times, aligned_spike_labels] = ...
    trialAlignSpikes(spike_times_cell, align_time, time_before, time_after); %FILLIN based on the function


%now we need to visualize this data by making a raster plot
%In lab0, you wrote a function plotRaster for making these plots. Use it
%here.

%loop through channels to make a figure + subplot with a raster for each
%channel. Save handles of figures so we can add subplot with the PSTH.
for iCh=1:num_channels
    fig_handle(iCh) = figure;
    ax_handle(iCh,1) = subplot(2,1,1);

    plotRaster(aligned_spike_times{iCh}, aligned_spike_labels{iCh}, fig_handle(iCh), ax_handle(iCh, 1)) %FiLLIN
    title(['Channel ', num2str(iCh) ' Stimulus response'])

end


%3. Next, we need to write code to calculate and plot a PSTH
% To do this, we will convert our spike-times into spike rates ('binning')
% Then, we can compute the average rate across trials (the PSTH).
%
%TECHNICAL NOTE: We could also 'bin' our spikes across all trials together
%(That is, we could do trial-averaging and binning simultaneously in a single calculation).
%These two ways of computing a PSTH are equivalent and equally useful for this dataset.
%In future data analysis (e.g. lab 4), we will explore data with more complex trial structures.
%In that case, binning first and then grouping trials for a PSTH gives a more compact data
%representation (spike rate matrix)that can be mined/explored more flexibly.
%
%We've provided a shell function binTrialAlignedSpikes.m for the binning
%step. Open it and fill it in, then use it here to calculate and plot a
%PSTH

bin_width = 0.01;
[trial_spike_rate, time_bins] = binTrialAlignedSpikes(aligned_spike_times, aligned_spike_labels, time_before, time_after, bin_width);

%get average across trials
trialAvg_spike_rate = mean(trial_spike_rate, 1); %FILLIN

%plot the PSHT along with the raster plots
for iCh=1:num_channels

    figure(fig_handle(iCh));
    ax_handle(iCh,2) = subplot(2,1,2);

    stairs(time_bins, trialAvg_spike_rate(:,:,iCh))

    xlabel('Time (s)')
    ylabel('Average Firing Rate')

end

%I recommend playing around with the bin size to understand how it
%influences our estimate of firing rate and the temporal dynamics.



%4. Finally, let's estimate our stimulus response and look at our
%trial-to-trial variability in response

%Compute the stimulus response by getting the firing rate change
%post-stimulus. Do this both for the trial-averaged response and on a
%single trial basis
baseline_window = [-1 0]; %FILLIN - pick a timewindow for estimating the baseline firing rate
baseline_idx = time_bins>=baseline_window(1) & time_bins<=baseline_window(2);

response_window = [0 1]; %FILLIN - pick a time window for estimating the response to tactile stimulus
response_idx = time_bins>=response_window(1) & time_bins<=response_window(2);

trAvg_response = squeeze(mean(trialAvg_spike_rate(:, response_idx, :), 2) - mean(trialAvg_spike_rate(:, baseline_idx, :), 2)); %FILLIN

singleTr_response = mean(trial_spike_rate(:, response_idx, :), 2) - mean(trial_spike_rate(:, baseline_idx, :), 2); %FILLIN

%plot a histogram of single trial responses to look at variabilityA
figure
hist(singleTr_response, 10)
hold on
plot([trAvg_response trAvg_response], [0 10], 'r--', 'lineWidth',2)
xlabel('Stimulus response')
ylabel('# trials')

%is your distribution reasonably gaussian?
%What are the main sources of variability from trial-to-trial?
%Why would we expect this variability to be gaussian if we have enough trials?




%%  Now we can proceed to look at trends in our data across our recordings

%start with a fresh workspace
clear all

%again, we want to point to where our data files are. Except now, we want to specify a list of all files we want to analyze.
dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 2\Lab 2 (Data)\'; %FILLIN with the path to where your data is stored

file_prefix = '2.'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type
events_file_type = '.txt';
% file_date   = ''; %FILLIN with the date string used in your file
file_idnums  = [1 2] ; %FILLIN with the LIST of numeric values of the file ids.
channels = [1 2];


%we also want to define the meta-data associated with each file we listed so we can analyze trends.
% stimulus_position = [1 2]; %FILLIN

threshold = [-0.08 -0.2; 0.06 0.1];
num_files = length(file_idnums);
num_channels = 2;

% loop through your files to generate PSTH, rasters, and stimulus response for each stimulus
% you presented.
for iF=1:num_files
    %FILLIN with relevant code (i.e. the same calculations we performed on the test file)
    full_file_name = [dataDir file_prefix num2str(file_idnums(iF)) file_type];
    events_file_name = [dataDir file_prefix num2str(file_idnums(iF)) events_file_type];
    
    [data, FS] = audioread(full_file_name);
    [EVENTS, EVENT_TIMES] = readEventsFile(events_file_name);

    spike_times_cell = cell(num_channels, 1);
    spike_waveforms_cell = cell(num_channels, 1);
    time_waveforms = [-10:19]./FS;

    for iCh=1:num_channels
        [spike_times_cell{iCh}, spike_waveforms_cell{iCh}] = detectSpikes(data(:, iCh), threshold(iF, iCh), FS);
        figure;
        plot(time_waveforms, spike_waveforms_cell{iCh}) %FILLIN with relevant variable name
        xlabel('Time (s)') %FILL IN label
        ylabel('Waveforms (mV)') %FILL IN label, be sure to include units
        set(gca, 'ylim') %keep axis the same across plots for easier comparison
        title(['Spike Waveform - Channel ', num2str(channels(iCh))])
        grid on

        % figure
        % time_data = [0:length(data(:,1))-1]./FS; %this creates a time vector of the data in seconds
        % plot(time_data, data(:,iCh)) %plot the raw data
        % xlabel('Time (s)')
        % ylabel('Voltage (mV)')
        % hold on
        % plot([spike_times_cell{iCh} spike_times_cell{iCh}]', repmat([.1 .15], length(spike_times_cell{iCh}),1)', 'r')

    end
    


    align_time = EVENT_TIMES;
    time_before = 1;
    time_after  = 1;
    [aligned_spike_times, aligned_spike_labels] = ...
        trialAlignSpikes(spike_times_cell, align_time, time_before, time_after); %FILLIN based on the function

    for iCh=1:num_channels
        fig_handle(iCh) = figure;
        ax_handle(iCh,1) = subplot(2,1,1);
    
        plotRaster(aligned_spike_times{iCh}, aligned_spike_labels{iCh}, fig_handle(iCh), ax_handle(iCh, 1)) %FiLLIN
        title(['Channel ', num2str(iCh) ' Stimulus response'])

    end
    


    bin_width = 0.01;
    [trial_spike_rate, time_bins] = binTrialAlignedSpikes(aligned_spike_times, aligned_spike_labels, time_before, time_after, bin_width);
    trialAvg_spike_rate = mean(trial_spike_rate, 1); %FILLIN

    for iCh=1:num_channels
        figure(fig_handle(iCh));
        ax_handle(iCh,2) = subplot(2,1,2);
    
        stairs(time_bins, trialAvg_spike_rate(:,:,iCh))
    
        xlabel('Time (s)')
        ylabel('Average Firing Rate')

    end


end %end loop through files.


% For your comprehension questions, attach figures of:
% Pick one example file and show the following raw data plots:
% - spike waveforms for channel 1
%  - spike waveforms for channel 2 (if you recorded 2 channels)
%  - trial-aligned raster and PSTH for channel 1
%  - trial-aligned raster and PSTH for channel 2 (if you recorded 2 channels)
