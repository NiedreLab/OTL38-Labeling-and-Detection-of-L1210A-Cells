function [non_coinc_peaks1, non_coinc_peaks2, coinc_pk_count, coinc_peaks1, coinc_peaks2] = removeCoincPeaks(peaks1, peaks2, time, nameValArgs)
%REMOVECOINCPEAKS Returns peaks structs with peaks coincident to eachother
% removed. Also cleans up each peak struct individually

arguments
    peaks1 struct
    peaks2 struct
    time(:,1) double
    
    nameValArgs.CoincidenceWindow double = 0.03;
end
coinc_window = nameValArgs.CoincidenceWindow;

fs = 1 ./ (time(2)-time(1)); % Sampling frequency
cws = coinc_window * fs;

% find any peaks that match each other between the two fibers within the
% user-specified coincidence window
matches1 = false(peaks1.count, 1);
matches2 = false(peaks2.count, 1);
for ii = 1:peaks1.count
    pk_loc = peaks1.locs(ii);
    coinc_loc = find(peaks2.locs>pk_loc-cws & peaks2.locs<pk_loc+cws);
    if(~isempty(coinc_loc))
        % if any are found indicate matches
        matches1(ii) = true;
        matches2(coinc_loc) = true;
    end
end

% remove all matches from peaks1
non_coinc_peaks1 = filterPeaks(peaks1, ~matches1);
non_coinc_peaks2 = filterPeaks(peaks2, ~matches2);
coinc_peaks1 = filterPeaks(peaks1, matches1);
coinc_peaks2 = filterPeaks(peaks2, matches2);

% Count the number of peaks removed
coinc_pk_count = sum(matches1) + sum(matches2);

end

