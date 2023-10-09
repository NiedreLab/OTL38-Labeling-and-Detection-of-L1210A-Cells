function [data_bs, noise, peaks, thresh_curve] = preProc1(data, time,...
    background_window, smooth_window, std_window, dynamic_window,...
    noise_alpha, std_thresh, rel_thresh, prominence_factor, plot_flag, params)
%PREPROC1 Pre-processes the data stored in the columns of the data argument
% Does background subtraction, estimates noise, calculates a moving peak 
% threshhold, and identifies peak candidates. 
% Can display information through different plots if enabled

% preProc1 emulates old fast background subtraction code, which is still
% kind of slow due to the way it handles noise calculation. 
% There's also a... confused implementation of the moving peak theshold 
% due to the limitations of the findpeaks(...) function

% 2023_01_31 - Amber - Added "round" to the fs calculation to avoid a
% rounding error in dynamic peak candidates when determining the number of
% intervals to use for the dynamic thresholding.

% Calculate basic values from user defined parameters
dt = time(2) - time(1);                     % time increment
fs = round(1 / dt);                                % sampling rate used in acquisition
bsw = background_window * fs;        % convert window to indicies
ssw = smooth_window * fs;                % convert window to indicies
stdsw = std_window * fs;                    % convert window to indicies
interval_len = dynamic_window * fs;         % convert window to indicies

% Smooth data to reduce noise
data_smooth = movmean(data, ssw);

% Basic median filter background subtraction
bg = movmedian(data_smooth, bsw);      % apply 'bsw' median filter,
data_bs =  data_smooth - bg;           % subtract background

% this sequence of code estimates the pre-processed signal background standard deviation (noise) over time.
std_sig = movstd(data_bs, stdsw);
std_proc = zeros(size(std_sig));
%std_int = false(size(std_smooth));

% Single Exponential Smoothing is used here because real peaks give a 
% transient increase in standard deviation, which should not be included 
% in the estimate of noise.
% Included/excluded regions are stored in std_int
std_proc(1:3,:) = std_sig(1:3,:);
for ii = 3:length(std_proc)
    % Exclude standard deviations as a result of peaks 
    for jj = 1:length(data(1,:))
        if std_sig(ii,jj) < std_thresh * std_proc(ii-1,jj)
            std_proc(ii,jj) = noise_alpha*std_sig(ii,jj) + (1-noise_alpha)*std_proc(ii-1,jj);
            %std_int(ii,jj) = true;
        else
            std_proc(ii,jj) = noise_alpha*std_proc(ii-1,jj) + (1-noise_alpha)*std_proc(ii-2,jj);
        end
    end
end

% Smooth std array again, for the hell of it
std_proc = movmean(std_proc, ssw);
noise = mean(std_proc);

% Generate a dynamic moving threshold from the processed movstd curve and
% find peak candidates using that threshold
[peaks, thresh_curve] = ...
    dynamicPeakCandidates(data_bs, std_proc, interval_len, rel_thresh, prominence_factor);

% Generate processing plots if so desired
if plot_flag
    preProcPlots(time, data, bg, data_bs, std_sig, std_proc, params)
end

end
