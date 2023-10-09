function [rm_bunched_peaks, bunched_pk_count] = removeBunchedPeaks(peaks, time, coinc_window)
%REMOVEBUNCHEDPEAKS Removes peaks
%   Detailed explanation goes here

cws = coinc_window ./ (time(2)-time(1)); % time window (s) * Hz (1/s)
rm_bunched_peaks = peaks;
bunched_pk_count = 0;

for ii = 1:length(peaks)
    good_pks = true(peaks(ii).count,1);
    for jj = 1:peaks(ii).count
        % If the current peak hasn't already been removed examine the
        % window around the current peak for bunched peaks
        if good_pks(jj)
            pk_loc = peaks(ii).locs(jj);
            % Find peaks coincident to eachother to within a coincidence window
            bunched_ind = find((peaks(ii).locs > pk_loc-cws) & (peaks(ii).locs < pk_loc+cws));
            % If there is more than one peak (that is, the current peak),...
            if length(bunched_ind) > 1
                % ... keep the largest peak and...
                [~, kk] = max(peaks(ii).pks(bunched_ind));
                bunched_ind = bunched_ind(bunched_ind ~= bunched_ind(kk));
                % ... remove the rest
                good_pks(bunched_ind) = false;
            end
        end
    end
    % remove all matches from peaks
    rm_bunched_peaks(ii) = filterPeaks(peaks(ii), good_pks);
    
    % Count the number of peaks removed
    bunched_pk_count = bunched_pk_count + sum(~good_pks);
end

end

