function [noise_flag1, noise_flag2] = checkNoisyPMT(data1, data2, time, params, nameValArgs)
% CHECKNOISYPMT Looks for obvious differences b/w 2 input sets of data from
% the same fiber.
% 
% Note that these are 2 PMTs from the same fiber, not 2 different fibers. 
% The shape of the curves should be very similar Since these are 
% measuring the same thing. Obvious variations are possilbly a sign of a 
% faulty PMT or system noise

% Validate arguments and parse name-value pairs
arguments
    data1 (:,1) double
    data2(:,1) double
    time (:,1) double
    params
    
    nameValArgs.SmoothWindow double = 0.01;
    nameValArgs.BackgroundWindow double = 5;
    nameValArgs.RelativeThresh double = 5;
    nameValArgs.ProminenceFactor double = 1;
    nameValArgs.NoiseThresh double = 3;
    nameValArgs.PeakColor (1,3) double = [1 0 1];
    nameValArgs.ThresholdColor (1,3) double = [0 1 1];
end

smooth_window = nameValArgs.SmoothWindow;
background_window = nameValArgs.BackgroundWindow;
rel_thresh = nameValArgs.RelativeThresh;
prominence_factor = nameValArgs.ProminenceFactor;
noise_thresh = nameValArgs.NoiseThresh;
pkColor = nameValArgs.PeakColor;
threshColor = nameValArgs.ThresholdColor;

% Relative threshold can be given as a scalar to be applied to all channels 
% or can be given as a vector of relative thresholds, each applied to a
% different source
if length(rel_thresh) == 1
    rel_thresh = ones(2,1) * rel_thresh;
elseif length(rel_thresh) ~= 2
    error('rel_thresh is neither a scalar nor is it equal in length to number of sources')
end

% Get the largest peaks for each data set 
[data_bs, ~, peaks, thresh_curve] = preProcLazy([data1 data2], time,...
    background_window, smooth_window,...
    rel_thresh, prominence_factor, false, []);

% A PMT is considered noisy if there are noise_thresh * # peaks more
% peaks in the current PMT than than the other PMT (If # peaks found
% in the other PMT is 0, look for noise_thresh # peaks in the current
% PMT
noise_flag1 = peaks(1).count > noise_thresh * peaks(2).count ||...
    (peaks(2).count == 0 && peaks(1).count > noise_thresh);
noise_flag2 = peaks(2).count > noise_thresh * peaks(1).count ||...
    (peaks(1).count == 0 && peaks(2).count > noise_thresh);

% Plot peaks from each channel for manual inspection
if ~isempty(params)
    peakThreshCurvePlot(data_bs, time, peaks, thresh_curve, params, 'o', ' PMT CHECK', pkColor, threshColor);
end

end

