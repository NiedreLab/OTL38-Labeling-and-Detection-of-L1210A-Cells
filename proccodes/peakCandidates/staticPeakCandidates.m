function [peaks, thresh_curve] = staticPeakCandidates(data_bs, peak_thresh, prominence_factor)
%STATICPEAKCANDIDATES 

sz = size(data_bs);

% peaks struct has 5 fields
% pks: Height of peaks
% locs: Location (index) of peaks
% widths: Widths of peaks
% proms: Prominence of each peak
% count: Number of peaks
peaks = struct('pks', [], 'locs', [], 'widths', [], 'proms', [],... 
    'count', num2cell(zeros(sz(2), 1)));

thresh_curve = ones(sz);
warning('off', 'signal:findpeaks:largeMinPeakHeight')
for ii = 1:sz(2)
    [peaks(ii).pks, peaks(ii).locs, peaks(ii).widths, peaks(ii).proms] = ...
        findpeaks(data_bs(:,ii),...
        'MinPeakHeight', peak_thresh(ii),...
        'MinPeakProminence', peak_thresh(ii) .* prominence_factor);
    peaks(ii).count = length(peaks(ii).pks);
    thresh_curve(:,ii) = peak_thresh(ii) .* thresh_curve(:,ii); 
end
warning('on', 'signal:findpeaks:largeMinPeakHeight')

end

