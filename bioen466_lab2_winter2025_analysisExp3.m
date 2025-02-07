%%%%% This script provides a guided outline for analyzing your experimental data collected for Experiment 3 (quantifying sensory response dynamics)
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

dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 2\Lab 2 (Data)\'; %FILLIN with the path to where your data is stored

file_prefix = '3.'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

% file_date   = ''; %FILLIN  with the date string used in your file
file_idnum  = [1 2]; %FILLIN with the numeric value of the file id

full_file_name = [dataDir file_prefix num2str(file_idnum) file_type];

% load your data file and the sampling rate into matlab into the variables 'data' and FS, respectively.
% hint: look at the Matlab function 'audioread'
[data, FS] = audioread(full_file_name);

% check the basic properties of your loaded variables
whos data FS

%note that if you are recording 2 channels instead of 1, your data will have 2 columns (#time points x #channels). Let's create a variable to keep track of how many channels we have so our code is flexible:
num_channels = %FILLIN

% plot 1 second of your data (channels 1 and 2) to assure the loaded file looks correct
figure
plot( data(1:FS,: ) ) %FILLIN the portion in the brackets to only plot 1 second (hint: you need to take the sampling rate into account)
xlabel('Time, in samples') %fill in the units of the time axis on this plot
ylabel('Voltage (mV)')


% Now let's find spikes and their waveforms in the data.
% use the code you developed in lab 1 for this task (detectSpikes.m)
%note that you may need to modify this function or how you call it if you have two channels of data in your file.
%create a variable spike_times - cell {# channels x 1} with the detected spike times for each channel.
%
%FILLIN with relevant code


%inspect your data to convince yourself the function is working as you expect
%Some of the plots we generated in lab 1 may be helpful here (e.g. raw
%traces + detected spike times)



%%%%%
%1. load an 'events' file.
events_file_type = '_events.txt';
events_file_name = [dataDir file_prefix file_date '_' num2str(file_idnum) events_file_type];

[EVENTS, EVENT_TIMES] = readEventsFile(events_file_name);

%IMPORTANT: double-check that you see two different types of events, one indicating
%stimulation start, and one indicating stimulation end.


%2. trial-align your spikes to the event corresponding to stimulation
%starting
align_event = ;%which code to look for in EVENTS
align_time = ; %use logical indexing to sub-select the EVENT_TIMES to use
time_before = ; %FILLIN
time_after  = ; %FILLIN

[aligned_spike_times, aligned_spike_labels] = ...
    trialAlignSpikes(); %FILLIN

%3. Compute the estimated rates
bin_width = 0.01;
[trial_spike_rate, time_bins] = binTrialAlignedSpikes(); %FILLIN

%get average across trials
trialAvg_spike_rate = ; %FILLIN

%loop through channels to make a figure + subplot with a raster for each
%channel and a PSTH.
for iCh=1:num_channels
    fig_handle = figure;
    ax_handle = subplot(2,1,1);

    plotRaster() %FILLIN
    title(['Channel ', num2str(iCh) ' Stimulus response'])



    subplot(2,1,2);
    stairs(time_bins, trialAvg_spike_rate(:,:,iCh))
    xlabel() %FILLIN
    ylabel() %FILLIN
end


%Now try repeating this procedure, but align to the time when you release
%the stimulation.




%%  Now we can proceed to look at trends in our data across our recordings

%start with a fresh workspace
clear all

%again, we want to point to where our data files are. Except now, we want to specify a list of all files we want to analyze.
dataDir = ''; %FILLIN  with the path to where your data is stored

file_prefix = 'BYB_Recording_'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

file_date   = ''; %FILLIN  with the date string used in your file
file_idnums  = [] ; %FILLIN with the LIST of numeric values of the file ids.


%we also want to define the meta-data associated with each file we listed so we can analyze trends.
stimulus_position = []; %FILLIN


num_files = length(file_idnums);

% loop through your files to generate PSTH, rasters, and stimulus response for each stimulus
% you presented.
for iF=1:num_files

  %FILLIN with relevant code (i.e. the same calculations we performed on the test file)

end %end loop through files.


% Finally, generate a figure to quantify and visualize the temporal dynamics of the
% sensory responses.
%Are response dynamics the same across neurons? Across positions?



% For your comprehension questions, attach figures of:
% Pick one example file and show the following raw data plots:
% - spike waveforms for channel 1
% - trial-aligned raster and PSTH for channel 1 or 2, aligned to stim on
% - trial-aligned raster and PSTH for channel 1 or 2, aligned to stim off
% - Summary figure(s) to visualize/quantify the temporal dynamics of the response (as described at end of script)
