function peakThreshCurvePlot(data, time, peaks, thresh_curve, params, marker, plotSuffix, pkColor, threshColor)

figure();
tiledlayout('flow');
for ii = 1:length(data(1,:))
    nexttile()
    hold on
    plot(time, data(:,ii), 'color', params(ii).color)
    plot(time(peaks(ii).locs), peaks(ii).pks, marker, 'MarkerSize', 6, 'Color', pkColor)
    plot(time, thresh_curve(:,ii), 'Color', threshColor)
    title([params(ii).name plotSuffix], 'Interpreter', 'none')
    xlabel('Time (s)')
    ylabel(['Signal (' params(ii).units ')'])
    legend('DiFC Data', 'Peaks', 'Threshold')
    hold off
end

end

