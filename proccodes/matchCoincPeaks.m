function [matched_peaks1, matched_peaks2, matched_score] =...
    matchCoincPeaks(peaks1, peaks2, time, nameValArgs)
%MATCHDIRECTIONALPEAKS Matches peaks in from the first set of peaks to
%peaks in the second channel

arguments
    peaks1 struct
    peaks2 struct
    time (:,1) double
    
    nameValArgs.CoincidenceWindow double = 0.03;
    nameValArgs.MaxScore double = 3;
end
coinc_window = nameValArgs.CoincidenceWindow; % seconds (s)
max_score = nameValArgs.MaxScore;  % Unitless

% Factor scaling heights in peaks2 during scoring. the ratio of the mean height of
% peaks1:peaks2 (peaks2 * peaks1/peaks2 ~= peaks1, where ~= is roughly
% equals LIKE IT SHOULD BE, MATLAB)
pks2Scale = mean(peaks1.pks) ./ mean(peaks2.pks);

% PeakTimes indicate when in the duration of the scan the peaks occured
peakTimes1 = time(peaks1.locs);
peakTimes2 = time(peaks2.locs);

% The match_struct is a list of all possible matches within the coincidence 
% and score thresholds. It stores the index of each peak from peaks1/2,
% respectively, and the score of the match
match_struct = repmat(struct('peak', 0, 'match', 0, 'score', 0), peaks1.count, 1);
kk = 1;

% Go through peaks1 and find all possible matches for each peak as
% determined by their location and score
for ii = 1:peaks1.count
    time_lo = peakTimes1(ii) - coinc_window;
    time_hi = peakTimes1(ii) + coinc_window;
    
    % Find peaks in the time range
    match_pk_ind = find(peakTimes2 > time_lo & peakTimes2 < time_hi);
    
    % If we find any match candidates, score them
    if ~isempty(match_pk_ind)
        score = scoreMatches(ii, match_pk_ind, peaks1, peaks2, pks2Scale);
        
        % Look for matches that score below the threshold. If there are
        % any, add the match to the match_struct.
        good_score = find(score < max_score)';
        for jj = good_score
            match_struct(kk).peak = ii;
            match_struct(kk).match = match_pk_ind(jj);
            match_struct(kk).score = score(jj);
            kk = kk + 1;
        end
    end
end

% Remove peaks with no possible matches / extra allocated structs
match_struct = match_struct([match_struct.peak] ~= 0);

% Sort match_struct by score, ascending
[~, idx] = sort([match_struct.score]);
match_struct = match_struct(idx);

% Create a bit map the length of the match_struct
best_matches = true(length(match_struct),1);
all_peaks = [match_struct.peak]';
all_matches = [match_struct.match]';

% Iterate throuhght the sorted match struct. We can declare matches with
% the lowest unmatched score the best, and use the bitmap to filter out any
% entries in match struct that contain either peak just declared to be the
% best
for ii = 1:length(match_struct)
    % If the current match has already been removed from contention move on
    if ~best_matches(ii); continue; end
    
    % Remove other entries in the match_struct contianing either peak
    best_matches = best_matches & ~(all_peaks == match_struct(ii).peak);
    best_matches = best_matches & ~(all_matches == match_struct(ii).match);
    
    % The removal process removes the current peak as well, so we set the
    % current peak to true
    best_matches(ii) = true;
end

% Filter the match_struct to keep the best matches
match_struct = match_struct(best_matches);

% Sort by the first peak
[~, idx] = sort([match_struct.peak]);
match_struct = match_struct(idx);

% Filter each array in both peaks structs by the best peaks/matches
all_best_peaks = [match_struct.peak];
all_best_matches = [match_struct.match];

% Return peaks structs with only matched peaks
matched_peaks1 = filterPeaks(peaks1, all_best_peaks);
matched_peaks2 = filterPeaks(peaks2, all_best_matches);

% Return matched speeds/scores
matched_score = [match_struct.score]';
end

% This function returns the scores for peak from peaks1 scored
% against all possible matches from peaks2
function score = scoreMatches(peak, matches, peaks1, peaks2, pks2Scale)

% Heights of peak and matches normalized to better resemble the height of
% the peak from peaks1
peakHeight = peaks1.pks(peak);
matchesHeight = peaks2.pks(matches) * pks2Scale;


% Widths of peak and matches in seconds (how long it took for each
% detection to pass the probe)
peakWidth = peaks1.widths(peak);
matchesWidth = peaks2.widths(matches);

% Here we score matched peaks (from peaks2) by how closely their
% widths and heights match the current peak in peaks1
% scoring ratios by abs(log2(ratio)) weights ratios of 1/2 as much
% as ratios of 2, both with a score of 1.
height_score = abs(log2(matchesHeight / peakHeight));
width_score = abs(log2(matchesWidth / peakWidth));

score = width_score + height_score;

end

