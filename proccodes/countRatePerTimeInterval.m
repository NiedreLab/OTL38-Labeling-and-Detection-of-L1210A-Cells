function [fwd_peaks_locs_in_time_interval, rev_peaks_locs_in_time_interval, peaks_P1_locs_in_time_interval,peaks_P2_locs_in_time_interval, time] = countRatePerTimeInterval(fwd_peaks,rev_peaks,peaks,time_interval, scan_length,sampling_rate)
%countRatePerTimeInterval-- Counts number of cells per specified time interval
%Josh Pace, 20230807
conversionFactor = sampling_rate*60; %To help convert peak clocations from sampleRate*seconds to minutes
time = 0:time_interval:(scan_length/60); %Setting up time based on the specified  time interval 


fwd_peaks_locs_in_min = fwd_peaks(1).locs/conversionFactor; % Locs in minutes 
fwd_peaks_locs_in_time_interval = [];

rev_peaks_locs_in_min = rev_peaks(1).locs/conversionFactor; % Locs in minutes 
rev_peaks_locs_in_time_interval = [];


peaks_P1_locs_in_min = peaks(1).locs/conversionFactor; % Locs in minutes 
peaks_P1_locs_in_time_interval = [];

peaks_P2_locs_in_min = peaks(2).locs/conversionFactor; % Locs in minutes 
peaks_P2_locs_in_time_interval = [];

for i = time
    fwd_peaks_locs_in_time_interval = [fwd_peaks_locs_in_time_interval numel(find(fwd_peaks_locs_in_min <= (i + time_interval) & fwd_peaks_locs_in_min >=i))];
    rev_peaks_locs_in_time_interval = [rev_peaks_locs_in_time_interval numel(find(rev_peaks_locs_in_min <= (i + time_interval) & rev_peaks_locs_in_min >=i))];
    peaks_P1_locs_in_time_interval = [peaks_P1_locs_in_time_interval numel(find(peaks_P1_locs_in_min <= (i + time_interval) & peaks_P1_locs_in_min >=i))];
    peaks_P2_locs_in_time_interval = [peaks_P2_locs_in_time_interval numel(find(peaks_P2_locs_in_min <= (i + time_interval) & peaks_P2_locs_in_min >=i))];
end
%Before a 0 was added at the start of the scan but after with discussing with mark, decided to stop that
fwd_peaks_locs_in_time_interval = [fwd_peaks_locs_in_time_interval(1) fwd_peaks_locs_in_time_interval(1:end-1)];
rev_peaks_locs_in_time_interval= [rev_peaks_locs_in_time_interval(1) rev_peaks_locs_in_time_interval(1:end-1)];
peaks_P1_locs_in_time_interval = [rev_peaks_locs_in_time_interval(1) peaks_P1_locs_in_time_interval(1:end-1)];
peaks_P2_locs_in_time_interval = [peaks_P2_locs_in_time_interval(1) peaks_P2_locs_in_time_interval(1:end-1)];

 
end