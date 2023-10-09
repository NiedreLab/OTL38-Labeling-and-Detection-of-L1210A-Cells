function [matched_peaks1, matched_peaks2, matched_speed, matched_score] =...
    matchDirectionalPeaks(peaks1, peaks2, time, probe_distance, nameValArgs)
%MATCHDIRECTIONALPEAKS Matches peaks in from the first set of peaks to
%peaks in the second channel

arguments
    peaks1 struct
    peaks2 struct
    time (:,1) double
    probe_distance double
    nameValArgs.MaxSpeed double = 300;
    nameValArgs.MaxSpeedFactor double = 10;
    nameValArgs.MaxScore double = 4;
    nameValArgs.Direction char = 'fwd';
end
max_speed = nameValArgs.MaxSpeed; % mm/s
min_speed_factor = 10; % Unitless
max_score = nameValArgs.MaxScore; % Unitless
direction = nameValArgs.Direction;

%  Only accepts fwd, rev as directions
if ~(strcmp(direction, 'fwd') || strcmp(direction, 'rev'))
    error('Invalid Direction')
end

% If direction is rev, swap input peaks
if strcmp(direction, 'rev')
    tmp = peaks2;
    peaks2 = peaks1;
    peaks1 = tmp;
end

% Calculate some useful values
dt = time(2) - time(1);
max_speed_time = probe_distance/max_speed; % mm / (mm/s) = s Min speed time in s

% Factor scaling heights in peaks2 during scoring. the ratio of the mean height of
% peaks1:peaks2 (peaks2 * peaks1/peaks2 ~= peaks1, where ~= is roughly
% equals LIKE IT SHOULD BE, MATLAB)
pks2Scale = mean(peaks1.pks) ./ mean(peaks2.pks);

% PeakTimes indicate when in the duration of the scan the peaks occured
peakTimes1 = time(peaks1.locs);
peakTimes2 = time(peaks2.locs);

% The match_struct is a list of all possible matches within the speed and
% score thresholds. It stores the index of each peak from peaks1/2,
% respectively, and the speed and score of the match
match_struct = repmat(struct('peak', 0, 'match', 0, 'score', 0, 'speed', 0), peaks1.count, 1);
kk = 1;

% Go through peaks1 and find all possible matches for each peak as
% determined by their location and score
for ii = 1:peaks1.count
    % We estimate the speed of the cell as 1/(the time it took the cell to
    % pass the probe), or 1 / (peak width)
    % est_speed = 1 ./ widths1(ii); % 1/ s
    % Use this estimation to estimate travel time
    % est_dt = probe_distance ./ est_speed; % mm / (1/s) = mm * s
    % in short:
    est_dt = probe_distance .* peaks1.widths(ii) * dt; % mm * s
    % Or estimating that the probe detects a 1mm long region in space,
    % mm / (mm/s) = s
    
    % Is this a reasonable time range?
    % Comment from the original function (I don't understand the it but ah well):
    % look for corresponding peak(s) in second channel from 0.1 s to up to 5X the estimated time delay
    time_lo = max_speed_time + peakTimes1(ii);
    time_hi = min_speed_factor * est_dt + peakTimes1(ii);
    
    % Find peaks in the time range
    match_pk_ind = find(peakTimes2 > time_lo & peakTimes2 < time_hi);
    
    % If we find any match candidates, score them
    if ~isempty(match_pk_ind)
        [score, speed] = scoreMatches(ii, match_pk_ind, peaks1, peaks2,...
            time, pks2Scale, dt, probe_distance);
        
        % Look for matches that score below the threshold. If there are
        % any, add the match to the match_struct.
        good_score = find(score < max_score)';
        for jj = good_score
            match_struct(kk).peak = ii;
            match_struct(kk).match = match_pk_ind(jj);
            match_struct(kk).score = score(jj);
            match_struct(kk).speed = speed(jj);
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

% If direction is rev, swap output peaks
if strcmp(direction, 'rev')
    tmp = matched_peaks2;
    matched_peaks2 = matched_peaks1;
    matched_peaks1 = tmp;
end

% Return matched speeds/scores
matched_speed = [match_struct.speed]';
matched_score = [match_struct.score]';

end

% This function returns the scores/speeds for peak from peaks1 scored
% against all possible matches from peaks2
function [score, speed_act] = scoreMatches(peak, matches, peaks1, peaks2, time, pks2Scale, dt, probe_distance)

% When in the duration of the scan the peaks occured
peakTime = time(peaks1.locs(peak));
matchesTime = time(peaks2.locs(matches));

% Heights of peak and matches normalized to better resemble the height of
% the peak from peaks1
peakHeight = peaks1.pks(peak);
matchesHeight = peaks2.pks(matches) * pks2Scale;


% Widths of peak and matches in seconds (how long it took for each
% detection to pass the probe)
peakWidth = peaks1.widths(peak) * dt;
matchesWidth = peaks2.widths(matches) * dt;

% Here we score matched peaks (from peaks2) by how closely their
% widths, heights and speeds match the current peak in peaks1
% scoring ratios by abs(log2(ratio)) weights ratios of 1/2 as much
% as ratios of 2, both with a score of 1.
height_score = abs(log2(matchesHeight / peakHeight));
width_score = abs(log2(matchesWidth / peakWidth));

% Find the actual speed of the matched peaks and the estimated
% speed of the matched peaks (as described above), and use their
% ratio as a score
speed_act = probe_distance ./ (matchesTime - peakTime); % mm / s
speed_est = 1 ./ ((matchesWidth+peakWidth)./2); % 1 / s (or mm / s)
speed_score = abs(log2(speed_act ./ speed_est));

score = width_score + height_score + speed_score;

end

