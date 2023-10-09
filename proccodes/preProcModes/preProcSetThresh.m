function [data_bs, noise, peaks, thresh_curve] = preProcSetThresh(data, time,...
    background_window, smooth_window, hard_thresh, prominence_factor,... 
    plot_flag, params)
%PREPROCSETTHRESH Pre-processes the data stored in the columns of the data argument
% Does background subtraction and identifies peak candidates.
% Can display information through different plots if enabled

% preProcSetThresh lets the user set their own threshhold for determining
% peak candidates

% Calculate basic values from user defined parameters
dt = time(2) - time(1);                     % time increment
fs = 1 / dt;                                % sampling rate used in acquisition
bsw = background_window * fs;       % convert window to datapoints
ssw = smooth_window * fs;         % convert window to datapoints

% Smooth data to reduce noise
data_smooth = movmean(data, ssw);

% Basic median filter background subtraction
bg = movmedian(data_smooth, bsw);      % apply 'bsw' median filter,
data_bs =  data_smooth - bg;                     % subtract background

% Estimate noise
noise = std(data_bs);

[peaks, thresh_curve] = staticPeakCandidates(data_bs, hard_thresh, prominence_factor);

if plot_flag
    preProcPlots(time, data, bg, data_bs, [], [], params)
end

