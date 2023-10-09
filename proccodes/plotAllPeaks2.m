function plotAllPeaks2(data, time, peaks, fwd_matches, rev_matches, thresh_curve, params)

threshColor = [180 58 58]./255;
unmatchedColor = [158 27 170]./255;
fwdColor = [38 137 187]./255;
revColor = [187 141 38]./255;

unmatched = peaks;
for ii = 1:2
    % build unmatched peak struct
    unmatched_arr = true(peaks(ii).count, 1);
    
    [~, jj, ~] = intersect(peaks(ii).locs, fwd_matches(ii).locs);
    unmatched_arr(jj) = false;
    
    [~, jj, ~] = intersect(peaks(ii).locs, rev_matches(ii).locs);
    unmatched_arr(jj) = false;
    
    unmatched(ii).pks = unmatched(ii).pks(unmatched_arr);
    unmatched(ii).locs = unmatched(ii).locs(unmatched_arr);
    unmatched(ii).widths = unmatched(ii).widths(unmatched_arr);
    unmatched(ii).proms = unmatched(ii).proms(unmatched_arr);
    
    % remove rev matches that overlap with fwd matches
    rev_arr = true(rev_matches(ii).count, 1);
    
    [~, jj, ~] = intersect(rev_matches(ii).locs, fwd_matches(ii).locs);
    rev_arr(jj) = false;
    
    rev_matches(ii).pks = rev_matches(ii).pks(rev_arr);
    rev_matches(ii).locs = rev_matches(ii).locs(rev_arr);
    rev_matches(ii).widths = rev_matches(ii).widths(rev_arr);
    rev_matches(ii).proms = rev_matches(ii).proms(rev_arr);
    
end

figure();
tiledlayout('flow');
for ii = 1:length(data(1,:))
    nexttile()
    hold on
    plot(time, data(:,ii), 'color', params(ii).color)
    plot(time, thresh_curve(:,ii), 'Color', threshColor)
    plot(time(unmatched(ii).locs), unmatched(ii).pks, 'o', 'MarkerSize', 6, 'Color', unmatchedColor)
    plot(time(fwd_matches(ii).locs), fwd_matches(ii).pks, '>', 'MarkerSize', 6, 'Color', fwdColor)
    plot(time(rev_matches(ii).locs), rev_matches(ii).pks, '<', 'MarkerSize', 6, 'Color', revColor)
    title([params(ii).name ' Peaks'], 'Interpreter', 'none')
    xlabel('Time (s)')
    ylabel(['Signal (' params(ii).units ')'])
    legend('DiFC Data', 'Threshold', 'Unmatched', 'Fwd', 'Rev')
    hold off
end

end

