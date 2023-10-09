function [peaks, thresh_curve] = dynamicPeakCandidates(data_bs, std_proc, interval_len, rel_thresh, prominence_factor)
%DYNAMICPEAKCANDIDATES

sz = size(data_bs);

% Generate a moving detection threshold
if mod(sz(1), interval_len) < 10^-5 % floating point hack
    scan_length_ints = floor(sz(1)/interval_len);  % The number of moving threshhold intervals in the scan
else
    scan_length_ints = ceil(sz(1)/interval_len);  % The number of moving threshhold intervals in the scan
end

lo = zeros(scan_length_ints, 1);
hi = zeros(scan_length_ints, 1);
moving_thresh = zeros(scan_length_ints, sz(2));
for ii = 1:scan_length_ints
    % lo and hi indicies of the current interval
    lo(ii) = int64((ii-1) * interval_len + 1);
    hi(ii) = int64(min(ii * interval_len, length(data_bs)));
    interval = std_proc(lo(ii):hi(ii),:);
    for jj = 1:sz(2)
        moving_thresh(ii,jj) = mean(interval(:,jj)) * rel_thresh(jj);
    end
end

% Smooth moving threshhold... once again for the hell of it...?
% movmean(A,5) = bsmooth(A)
moving_thresh = movmean(movmean(moving_thresh, 5, 1), 5, 1);

% peaks struct has 5 fields
% pks: Height of peaks
% locs: Location (index) of peaks
% widths: Widths of peaks
% proms: Prominence of each peak
% count: Number of peaks
peaks = struct('pks', [], 'locs', [], 'widths', [], 'proms', [],...
    'count', num2cell(zeros(sz(2), 1)));

% search for peaks using the adaptive threshold in each interval
thresh_curve = ones(size(data_bs));

warning('off', 'signal:findpeaks:largeMinPeakHeight')
for ii = 1:scan_length_ints
    for jj = 1:sz(2)
        peak_thresh = moving_thresh(ii, jj);
        [pks, locs, widths, proms] = findpeaks(data_bs(lo(ii):hi(ii),jj),...
            'MinPeakHeight', peak_thresh,...
            'MinPeakProminence', peak_thresh .* prominence_factor);
        % Update pks and locs data
        peaks(jj).pks = [peaks(jj).pks; pks];
        peaks(jj).locs = [peaks(jj).locs; locs+double(lo(ii))-1];
        peaks(jj).widths = [peaks(jj).widths; widths];
        peaks(jj).proms = [peaks(jj).proms; proms];
        peaks(jj).count = peaks(jj).count + length(pks);
    end
    thresh_curve(lo(ii):hi(ii),:) = thresh_curve(lo(ii):hi(ii),:) .* moving_thresh(ii,:);
end
warning('on', 'signal:findpeaks:largeMinPeakHeight')

end

