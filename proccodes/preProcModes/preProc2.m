function [data_bs, noise, peaks, thresh_curve] = preProc2(data, time,...
    background_window, smooth_window, std_window, dynamic_window,...
    noise_window, rel_thresh, prominence_factor, plot_flag, params)
%PREPROC2 Pre-processes the data stored in the columns of the data argument
% Does background subtraction, estimates noise, calculates a moving peak 
% threshhold, and identifies peak candidates. 
% Can display information through different plots if enabled

% preProc2 uses a median filter to exclude parts of the movstd curve that 
% are the result of noise. Uses the same confused implementation of the
% moving peak theshold arsing from the limitations of the 
% findpeaks(...) function

% Calculate basic values from user defined parameters
dt = time(2) - time(1);                     % time increment
fs = 1 / dt;                                % sampling rate used in acquisition
bsw = background_window * fs;       % convert window to datapoints
ssw = smooth_window * fs;                  % convert window to datapoints
snrsw = std_window * fs;                   % convert window to datapoints
nsw = noise_window * fs;
interval_len = dynamic_window * fs;           % convert window to datapoints

% Smooth data to reduce noise
data_smooth = movmean(data, ssw);

% Basic median filter background subtraction
bg = movmedian(data_smooth, bsw);              % apply 'bsw' median filter,
data_bs =  data_smooth - bg;                   % subtract background

% this sequence of code estimates the pre-processed signal background 
% standard deviation (noise) over time.
std_sig = movstd(data_bs, snrsw);

% Smooth moving std curve
std_smooth = movmean(std_sig, ssw);

% Apply the median filter to try to exclude std values resulting from noise
std_med = movmedian(std_smooth, nsw);

% Smooth std array again, for the hell of it and calculate noise
std_med = movmean(std_med, ssw);
noise = mean(std_med);

% Generate a dynamic moving threshold from the processed movstd curve and
% find peak candidates using that threshold
[peaks, thresh_curve] = ...
    dynamicPeakCandidates(data_bs, std_med, interval_len, rel_thresh, prominence_factor);

% Generate processing plots if so desired
if plot_flag
    preProcPlots(time, data, bg, data_bs, std_sig, std_med, params)
end

end
