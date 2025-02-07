dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 2\Lab 2 (Data)\'; %FILLIN in with the path to where your data is stored

file_prefix = '3.'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

% file_date   = ''; %FILLIN with the date string used in your file
file_idnum  = 2; %FILLIN with the numeric value of the file id

full_file_name = [dataDir file_prefix num2str(file_idnum) file_type];


% Load the .wav file
[data, FS] = audioread(full_file_name);  % Replace '3.2.wav' with your file path

% Use only the first column of the data
data = data(:, 1);

% Compute mean and standard deviation
mu = mean(data); % Mean of the data
sigma = std(data); % Standard deviation of the data

% Detect spikes using the program we wrote for experiment 1.
% Define a threshold as mean ± k * sigma
k = 2; % Multiplier (commonly 2 or 3)
threshold_upper = mu + k * sigma; % Upper threshold
threshold_lower = mu - k * sigma; % Lower threshold

threshold = threshold_upper; %FILLIN with a number based on your data.
[spike_times, spike_waveforms] = detectSpikes(data(:, 1), threshold(1), FS);

% Initialize the spike rate array
spike_rate = [];

% Loop through the data and calculate spike rate for each window
for i = 0 : 1 : length(spike_train)/FS - 1
    spike_idx = (spike_times >= i & spike_times <= i + 1); %FILLIN
    spike_count = sum(spike_idx);  % Count the spikes in the window
    spike_rate = [spike_rate; spike_count];  % Normalize by window size (Hz)
end

% Create the time axis for plotting
time_axis = (0:length(spike_rate)-1)';  % Each point is 1 second

% Plot the Time vs. Spike Rate
figure;
plot(time_axis, spike_rate, 'b');
xlabel('Time (ms)');
ylabel('Spike Rate (Hz)');
title('Time vs. Spike Rate');
grid on;