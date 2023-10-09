function [noise, newBackground] = getBackgroundAndNoise(data,multiplier)
%getBackgroundNoise: Calculates the noise based on the background
%   Background can be raised by 'multiplier' to simulate a higher signal and see the resulting noise

bsw = 5 * 2000; %background_window * fs;        % convert window to indicies
ssw = 0.01 * 2000;%smooth_window * fs;                % convert window to indicies
stdsw = 1 * 2000; %std_window * fs;                    % convert window to indicies
std_thresh = 1.05;
noise_alpha = 0.5;
if data < 0
    data = -data;
else
    data = data;
end

data = data*multiplier;%Simulate a higher background (i.e. OTL38 in circulation)

newBackground = data;
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
noise = mean(std_proc);
end