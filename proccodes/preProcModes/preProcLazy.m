function [data_bs, noise, peaks, thresh_curve] = preProcLazy(data, time,...
    background_window, smooth_window, rel_thresh, prominence_factor,... 
    plot_flag, params)
%PREPROCLAZY Pre-processes the data stored in the columns of the data argument
% Does background subtraction, estimates noise, and identifies peak candidates.
% Can display information through different plots if enabled

% preProcLazy uses the standard deviation to estimate the noise which saves
% loads of time spent calculating the moving threshold

% Calculate basic values from user defined parameters
sz = size(data);                            % size of input data
dt = time(2) - time(1);                     % time increment
fs = 1 / dt;                                % sampling rate used in acquisition
bsw = background_window * fs;               % convert window to indicies
ssw = smooth_window * fs;                % convert window to indicies

% Smooth data to reduce noise
data_smooth = movmean(data, ssw);

% Basic median filter background subtraction
bg = movmedian(data_smooth, bsw);               % apply 'bsw' median filter
data_bs =  data_smooth - bg;                    % subtract background

% Estimate noise
noise = std(data_bs);

% Peak threshhold is the estimated noise * declared relative threshhold
for ii = 1:sz(2)
    peak_thresh = noise .* rel_thresh(ii);
end

[peaks, thresh_curve] = staticPeakCandidates(data_bs, peak_thresh, prominence_factor);

if plot_flag
    preProcPlots(time, data, bg, data_bs, [], [], params)
end

