function [fig_handle, ax_handle] = plotRaster(aligned_spike_times, aligned_spike_labels, fig_handle, ax_handle)

%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
% The function description below described the high-level goals of the function and formats of the inputs and outputs. Read this carefully.

%[fig_handle] = plotRaster(aligned_spike_times, aligned_spike_labels)
%   Function to make a raster plot of trial-aligned spike times
%   Only plots a figure for a SINGLE channel
%   If no figure/axis provided, creates a new figure and axis.
%   If a figure and/or axis is provided, will draw raster there.
%
%inputs: aligned_spike_times - vector [#spikes x 1] with spike time-stamps, trial-aligned
%        aligned_spike_labels - vector [#spikes x 1] with labels (trial #) for each spike
%        fig_handle - handle to the figure to use for plot (optional. If
%                     not provided or empty, will make a new figure and axes)
%        ax_handle - handle to the axes to use for plot (optional. If not
%                    provided or empty, will make new axes)
%outputs: fig_handle - handle to the figure created with the raster plot.
%         ax_handle  - handle to the axes created with the raster plot.

%constants
TICK_HEIGHT = 0.5;

%figure handle logistics.
if ~exist('fig_handle', 'var') || isempty(fig_handle) %if no figure/axes, make one
    fig_handle = figure;
    ax_handle = axes;
elseif ~exist('ax_handle', 'var') || isempty(ax_handle) %if figure but no axes, draw them
    ax_handle = axes(fig_handle);
end
set(fig_handle, 'CurrentAxes', ax_handle) %set fig_handle, ax_handle to current axes. Plots will get drawn there.
hold on


%figure out how many trials exist
trial_list = unique(aligned_spike_labels);
num_trials = length(trial_list);

for iT=trial_list'

    %find spikes that occur on this trial
    %hint: this is what aligned_spike_labels tracks.
    idx = aligned_spike_labels == iT; %logical vector. %FILL IN
    num_spikes = sum(idx);

    %create X values for each tick mark
    %should be a [2 x #spikes] vector
    X = [aligned_spike_times(idx)'; aligned_spike_times(idx)']; %FILLIN

    %create Y values for each tick mark
    %shold be a [2 x #spikes] vector.
    %Y value for a spike should be [trial# trial#+TICK_HEIGHT]
    Y = [ones(1, num_spikes)*iT; ones(1, num_spikes)*iT+TICK_HEIGHT]; %FILLIN

    %plot your [X Y] lines
    plot(X,Y, 'r-')

end

%also plot a vertical line at X=0 to visualize alignment
plot([0 0], [0 num_trials+TICK_HEIGHT], 'k--', 'linewidth', 2)

%label axes
xlabel('Time (s)') %FILL IN
ylabel('Trial #') %FILL IN
