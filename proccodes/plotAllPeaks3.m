function plotAllPeaks3(data, time, peaks, fwd_matches, rev_matches,coinc_peaks, thresh_curve, params,saveDir,phantom,yAxesMax,yAxesMin)
%Adapted Malcolm's 'plotAllPeaks2'
%Modified by Josh Pace
%08/06/2021


warning('off','MATLAB:hgconvertunits:InvalidRefObj');

axisFontSize = 25;
labelFontSize = 30;
legendLineWidth = 2;
threshColor = [0 0 0]./255;
unmatchedColor = [158 27 170]./255;
fwdColor = [34 139 34]./255;
revColor = [238 0 0]./255;
coincColor = [1 1 0]./255;
%Determing if processing on a windows or mac
if ispc
    slash = '\';
else
    slash = '/';
end

if params(1).units == 'mV'
        systemName = 'b';
%     else
%         systemName = 'Red';
%         params(1).color = 'r';
%  end
else
        systemName = 'NIR';
         params(1).color = 'r';
 end



    sources = length(data(1,:));
unmatched = peaks;
for ii = 1:sources
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
    
%     % remove rev matches that overlap with fwd matches
%     rev_arr = true(rev_matches(ii).count, 1);
%     
%     [~, jj, ~] = intersect(rev_matches(ii).locs, fwd_matches(ii).locs);
%     rev_arr(jj) = false;
%     
%     rev_matches(ii).pks = rev_matches(ii).pks(rev_arr);
%     rev_matches(ii).locs = rev_matches(ii).locs(rev_arr);
%     rev_matches(ii).widths = rev_matches(ii).widths(rev_arr);
%     rev_matches(ii).proms = rev_matches(ii).proms(rev_arr);
%     
end
if sources == 2
matchedTotal = fwd_matches(1).count + rev_matches(1).count;

unMatchedF1 = peaks(1).count-matchedTotal;
if unMatchedF1 < 0
    unMatchedF1 = 0;
end
unMatchedF2 = peaks(2).count-matchedTotal;

if unMatchedF2 < 0
    unMatchedF2 = 0;
end
unmatchedTotal = [unMatchedF1 unMatchedF2];


if sources == 4
matchedTotalGFP = fwd_matches(1).count + rev_matches(1).count;
matchedTotaltdTomato = fwd_matches(2).count + rev_matches(2).count;

unMatchedP1GFP = peaks(1).count-matchedTotalGFP;
if unMatchedP1GFP < 0
    unMatchedP1GFP = 0;
end

unMatchedP2GFP = peaks(3).count-matchedTotalGFP;
if unMatchedP2GFP < 0
    unMatchedP2GFP = 0;
end
unMatchedP1tdTomato = peaks(2).count-matchedTotaltdTomato;

if unMatchedP1tdTomato < 0
    unMatchedP1tdTomato = 0;
end

unMatchedP2tdTomato = peaks(4).count-matchedTotaltdTomato;

if unMatchedP2tdTomato < 0
    unMatchedP2tdTomato = 0;
end
unmatchedTotal = [unMatchedP1GFP unMatchedP1tdTomato unMatchedP2GFP unMatchedP2tdTomato];
end

% Getting y axis limts so both probes have the same y axis
% dataMin = min(min(data));
% dataMax = max(max(data));
% yAxesMaxTmp = round(dataMax,-2);
% 
% if yAxesMaxTmp == 0
%         yAxesMax = 50;
%     else
%         if yAxesMaxTmp < dataMax
%             yAxesMax = 100 + yAxesMaxTmp;
%         else
%             yAxesMax = yAxesMaxTmp;
%         end
% end
% if le(dataMin, -50)
%     yAxesMin = -100;
% else
%     yAxesMin = -50;
% end
% 
% %     yAxesMax = 200;
% %     yAxesMin = -200;


if phantom == 1%For plotting phantom data
    figure(9);
    tiledlayout(sources,1);
    for ii = 1:length(data(1,:))
        nexttile()
        plot(time, data(:,ii), 'color', params(ii).color)
        hold on
        plot(time, thresh_curve(:,ii), 'Color', threshColor)
        plot(time(unmatched(ii).locs), unmatched(ii).pks, 'o', 'MarkerSize', 6, 'Color', unmatchedColor)
        plot(time(fwd_matches(ii).locs), fwd_matches(ii).pks, '>', 'MarkerSize', 6, 'Color', fwdColor)
        plot(time(rev_matches(ii).locs), rev_matches(ii).pks, '<', 'MarkerSize', 6, 'Color', revColor)
        title([params(ii).name], 'Interpreter', 'none','FontSize', 18)
        %ylabel(['Signal (' params(ii).units ')'])
        hold off
        legend('DiFC Data', 'Threshold', sprintf('Unmatched: %g',unmatchedTotal(ii)),sprintf('Fwd: %g',fwd_matches(ii).count),sprintf('Rev: %g',rev_matches(ii).count),'FontSize', 18)
        ylim([yAxesMin yAxesMax])
        axis = gca;
        axis.FontSize = axisFontSize;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        %ylabel( [systemName(ii) '-DiFC Signal (' params(ii).units ')'], 'FontSize',18,'FontWeight','bold')
    end
     fig = gcf;
     tl= fig.Children;
     xlabel(tl,'Time (s)','FontSize',labelFontSize,'FontWeight','bold');
     ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',labelFontSize,'FontWeight','bold')
     saveas(figure(9),[saveDir slash 'All Peaks.fig']); 
else
    figure(9);
    tiledlayout(sources,1);
    for ii = 1:length(data(1,:))
        nexttile()
        plot(time, data(:,ii), 'color', params(ii).color)
        hold on
        plot(time, thresh_curve(:,ii), 'Color', threshColor)
        plot(time(unmatched(ii).locs), unmatched(ii).pks, 'o', 'MarkerSize', 6, 'Color', unmatchedColor)
        plot(time(fwd_matches(ii).locs), fwd_matches(ii).pks, '>', 'MarkerSize', 6, 'Color', fwdColor)
        plot(time(rev_matches(ii).locs), rev_matches(ii).pks, '<', 'MarkerSize', 6, 'Color', revColor)
        plot(time(coinc_peaks(ii).locs), coinc_peaks(ii).pks, 'x', 'MarkerSize', 6, 'Color', coincColor)
        title([params(ii).name], 'Interpreter', 'none','FontSize', 18)
        %ylabel(['Signal (' params(ii).units ')'])
        %xlabel('Time (s)')
        legend('DiFC Data', 'Threshold', sprintf('Unmatched: %g',unmatchedTotal(ii)),sprintf('Fwd: %g',fwd_matches(ii).count), sprintf('Rev: %g',rev_matches(ii).count), sprintf('Coincident: %g',coinc_peaks(ii).count),'FontSize', 18,'Location','northeastoutside');
        axis = gca;
        axis.FontSize = axisFontSize;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        hold off
        ylim([yAxesMin yAxesMax])
        
    end
     fig = gcf;
     tl= fig.Children;
     xlabel(tl,'Time (s)','FontSize',labelFontSize,'FontWeight','bold');
     ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',labelFontSize,'FontWeight','bold')
     saveas(figure(9),[saveDir slash 'All Peaks.fig']); 
end
end
