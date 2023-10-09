function [data_bs, noise, peaks, thresh_curve, in_dat] = preProc(data, time, params, nameValArgs)
%PREPROC Wrapper Function for all the Preprocessing functions

% Validate arguments and parse name-value pairs
arguments
    data (:,:) double
    time (:,1) double
    params (:,1) struct = [];
    
    nameValArgs.Mode string = 'stdSmooth';
    nameValArgs.StartTrim double = 0;
    nameValArgs.EndTrim double = 0;
    nameValArgs.SmoothWindow double = 0.01;
    nameValArgs.BackgroundWindow double = 5;
    nameValArgs.HardThresh (:,1) double = nan;
    nameValArgs.RelativeThresh (:,1) double = 5;
    nameValArgs.ProminenceFactor double = 1;
    nameValArgs.StdWindow double = 1;
    nameValArgs.DynamicWindow double = 60;
    nameValArgs.NoiseAlpha double = 0.5;
    nameValArgs.StdThresh double = 1.05;
    nameValArgs.NoiseWindow double = 10;
    nameValArgs.RemoveBunchedPeaks logical = true;
    nameValArgs.CoincidenceWindow double = 0.03;
    nameValArgs.ShowDetailedPlots logical = false;
    nameValArgs.PeakColor (1,3) double = [1 0 1];
    nameValArgs.ThresholdColor (1,3) double = [0 1 1];
    %*********************************************************************%
    % If new name-value pairs are added, 
    % nameValArgs.NewName (Size) Class = Value;
    % goes here
    %*********************************************************************%
end

mode = nameValArgs.Mode;
start_trim = nameValArgs.StartTrim;
end_trim = nameValArgs.EndTrim;
smooth_window = nameValArgs.SmoothWindow;
background_window = nameValArgs.BackgroundWindow;
hard_thresh = nameValArgs.HardThresh;
rel_thresh = nameValArgs.RelativeThresh;
prominence_factor = nameValArgs.ProminenceFactor;
std_window = nameValArgs.StdWindow;
dynamic_window = nameValArgs.DynamicWindow;
noise_alpha = nameValArgs.NoiseAlpha;
std_thresh = nameValArgs.StdThresh;
noise_window = nameValArgs.NoiseWindow;
rm_bunched_flag = nameValArgs.RemoveBunchedPeaks;
coinc_window = nameValArgs.CoincidenceWindow;
plot_flag = nameValArgs.ShowDetailedPlots;
pkColor = nameValArgs.PeakColor;
threshColor = nameValArgs.ThresholdColor;
%*************************************************************************%
% If new name-value pairs are added, 
% value = nameValArgs.NewName;
% goes here
%*************************************************************************%

% Trim data
dt = time(2) - time(1);
fs = 1/dt;
lo = start_trim * fs;
hi = end_trim * fs;
time = time(lo+1:end-hi);
time = time - time(1);
data = data(lo+1:end-hi,:);
sz = size(data);


% If the data is on average negative, with negative peaks, it needs to be
% negated before pre-processing
if mean(data) < 0
    data = -data;
end

% Relative threshold can be given as a scalar to be applied to all channels 
% or can be given as a vector of relative thresholds, each applied to a
% different source
if length(rel_thresh) == 1
    rel_thresh = ones(sz(2),1) * rel_thresh;
elseif length(rel_thresh) ~= sz(2)
    error('rel_thresh is neither a scalar nor is it equal in length to number of sources')
end

% While we're at it, check for a valid hard_thresh as well...
if ~any(isnan(hard_thresh)) && (length(hard_thresh) ~= sz(2))
     error('hard_thresh is not equal in length to number of sources')
end

% If for some reason the user requests detailed plots but does not give
% a params struct, throw a warning and proceed without showing detailed
% plots
if isempty(params) && plot_flag
    warning('Detailed plots requested but given empty params struct.\n Detailed Plots will not be displayed')
    plot_flag = false;
end

% Turn off warning for too high of a threshold
warning('off', 'signal:findpeaks:largeMinPeakHeight')

%-------------------------------------------------------------------------%
% Choose Pre Processing Mode and Pre-process Data
%-------------------------------------------------------------------------%

if strcmp(mode, 'stdSmooth')
    % Pre-process using stdSmooth mode
    [data_bs, noise, peaks, thresh_curve] = preProc1(data, time,...
        background_window, smooth_window, std_window, dynamic_window,...
        noise_alpha, std_thresh, rel_thresh, prominence_factor, plot_flag, params);
    
    % Create in_dat struct for stdSmooth mode
    in_dat.mode = mode;
    in_dat.background_window = background_window;
    in_dat.smooth_window = smooth_window;
    in_dat.std_window = std_window;
    in_dat.dynamic_window = dynamic_window;
    in_dat.noise_alpha = noise_alpha;
    in_dat.std_thresh = std_thresh;
    in_dat.rel_thresh = rel_thresh;
    in_dat.prominence_factor = prominence_factor;
    
elseif strcmp(mode, 'stdMed')
    
    % Pre-process using stdMed mode
    [data_bs, noise, peaks, thresh_curve] = preProc2(data, time,...
        background_window, smooth_window, std_window, dynamic_window,...
        noise_window, rel_thresh, prominence_factor, plot_flag, params);
    
    % Create in_dat struct for stdMed mode
    in_dat.mode = mode;
    in_dat.background_window = background_window;
    in_dat.smooth_window = smooth_window;
    in_dat.std_window = std_window;
    in_dat.dynamic_window = dynamic_window;
    in_dat.noise_window = noise_window;
    in_dat.rel_thresh = rel_thresh;
    in_dat.prominence_factor = prominence_factor;
    
elseif strcmp(mode, 'lazyThresh')
    
    % Pre-process using lazyThresh mode
    [data_bs, noise, peaks, thresh_curve] = preProcLazy(data, time,...
        background_window, smooth_window, rel_thresh, prominence_factor,...
        plot_flag, params);
    
    % Create in_dat struct for lazyThresh mode
    in_dat.mode = mode;
    in_dat.background_window = background_window;
    in_dat.smooth_window = smooth_window;
    in_dat.rel_thresh = rel_thresh;
    in_dat.prominence_factor = prominence_factor;
    
elseif strcmp(mode, 'setThresh')
    
    % If the mode is setThresh, the HardThresh Name-Value pair MUST exist
    % That is, there is no default value
    if isnan(hard_thresh)
        error('Must set hard thresholds for setThresh Mode');
    end
    
    % Pre-process using setThresh mode
    [data_bs, noise, peaks, thresh_curve] = preProcSetThresh(data, time,...
        background_window, smooth_window, hard_thresh, prominence_factor,...
        plot_flag, params);
    
    % Create in_dat struct for setThresh mode
    in_dat.mode = mode;
    in_dat.background_window = background_window;
    in_dat.smooth_window = smooth_window;
    in_dat.hard_thresh = hard_thresh;
    in_dat.prominence_factor = prominence_factor;
    
    %*********************************************************************%
    % If new preProc modes are added, 
    % elseif strcmp(mode, 'newMode') 
    % ...
    % goes here
    %*********************************************************************%
    
else
    % If the mode is none of the above, reset threshold warning and 
    % throw an error
    warning('on', 'signal:findpeaks:largeMinPeakHeight')
    error('Invalid pre-processing mode')
end

% Reset threshold warning
warning('on', 'signal:findpeaks:largeMinPeakHeight')

% Remove bunched peaks if so desired
if rm_bunched_flag

    [peaks, ~] = removeBunchedPeaks(peaks, time, coinc_window);
    
    % Add the coincidence window to in_dat
    in_dat.coinc_window = coinc_window;
end

% Plot peak candidates if given params
if ~isempty(params)
    peakThreshCurvePlot(data_bs, time, peaks, thresh_curve, params, 'o', ' Peak Candidates', pkColor, threshColor)
end
end

