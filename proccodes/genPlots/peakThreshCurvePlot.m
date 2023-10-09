function peakThreshCurvePlot(data, time, peaks, thresh_curve, params, marker, plotSuffix, pkColor, threshColor)

figure();
tiledlayout(2,1);
for ii = 1:length(data(1,:))
    nexttile()
    hold on
    plot(time, data(:,ii), 'color', params(ii).color)
    plot(time(peaks(ii).locs), peaks(ii).pks, marker, 'MarkerSize', 6, 'Color', pkColor)
    plot(time, thresh_curve(:,ii), 'Color', threshColor)
    title([params(ii).name plotSuffix], 'Interpreter', 'none')
    xlabel('Time (s)')
    ylabel(['Signal (' params(ii).units ')'])
    legend('DiFC Data', 'Peaks', 'Threshold', 'Location', 'bestoutside')
    hold off
end

end

