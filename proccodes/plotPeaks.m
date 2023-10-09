function plotPeaks(data, time, peaks, thresh_curve, params, nameValArgs)
%PLOTPEAKS Summary of this function goes here
%   Detailed explanation goes here

arguments
    data (:,:) double
    time (:,1) double
    peaks(:,1) struct
    thresh_curve(:,:) double
    params (:,1) struct
    nameValArgs.Direction string = 'none';
    nameValArgs.PeakColor (1,3) double = [1 0 1];
    nameValArgs.ThresholdColor (1,3) double = [0 1 1];
end
direction = nameValArgs.Direction;
pkColor = nameValArgs.PeakColor;
threshColor = nameValArgs.ThresholdColor;

% Get the marker type for the match plot
if strcmp(direction, 'fwd')
    marker = '>';
    plotSuffix = ' Forward Matched Peaks'; 
elseif strcmp(direction, 'rev')
    marker = '<';
    plotSuffix = ' Reverse Matched Peaks'; 
elseif strcmp(direction, 'coinc')
    marker = 'x';
    plotSuffix = ' Coincident Peaks'; 
elseif strcmp(direction, 'none')
    marker = 'o';
    plotSuffix = ' Peaks'; 
else
    error('Invalid Direction')
end

peakThreshCurvePlot(data, time, peaks, thresh_curve, params, marker, plotSuffix, pkColor, threshColor);
end

