dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 2\Lab 2 (Data)\'; %FILLIN in with the path to where your data is stored

file_prefix = '3.'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

% file_date   = ''; %FILLIN with the date string used in your file
file_idnum  = 2; %FILLIN with the numeric value of the file id

full_file_name = [dataDir file_prefix num2str(file_idnum) file_type];

% Load the .wav file
[spike_train, fs] = audioread(full_file_name);  % Replace '3.2.wav' with your file path

% Use only the first column of the data
spike_train = spike_train(:, 1);

% Define the window size (1 second) and step size (1 second)
window_size = fs;  % 1 second window (fs is the sampling frequency)
step_size = fs;    % Move the window by 1 second

% Initialize the spike rate array
spike_rate = [];

% Loop through the data and calculate spike rate for each window
for i = 1:step_size:length(spike_train)-window_size
    window = spike_train(i:i+window_size-1);
    spike_count = sum(window);  % Count the spikes in the window
    spike_rate = [spike_rate; spike_count / window_size];  % Normalize by window size (Hz)
end

% Create the time axis for plotting
time_axis = (0:length(spike_rate)-1)';  % Each point is 1 second

% Plot the Time vs. Spike Rate
figure;
plot(time_axis, spike_rate, 'b');
xlabel('Time (s)');
ylabel('Spike Rate (Hz)');
title('Time vs. Spike Rate');
grid on;