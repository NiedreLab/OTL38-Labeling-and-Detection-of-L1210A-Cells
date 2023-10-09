function seePeaks(data, time, peaks, thresh_curve, params, nameValArgs)

arguments
    data (:,:) double
    time (:,1) double
    peaks(:,1) struct
    thresh_curve(:,:) double
    params (:,1) struct
    
    nameValArgs.Channel double = 1;
    nameValArgs.StartTime = NaN;
    nameValArgs.StartLoc = NaN;
    nameValArgs.PeakColor (1,3) double = [1 0 0];
    nameValArgs.ThresholdColor (1,3) double = [0 0 0];
    nameValArgs.MarkerColor (1,3) double = [.380 0 0];
end

channel = nameValArgs.Channel;
startTime = nameValArgs.StartTime;
startLoc = nameValArgs.StartLoc;
peakColor = nameValArgs.PeakColor;
threshColor = nameValArgs.ThresholdColor;
markerColor = nameValArgs.MarkerColor;


sz = size(data);
if channel > sz(2)
    error('Channel %g is larger than number of channels in data %g', channel, sz(2))
end

if ~isnan(startTime) && ~isnan(startLoc)
    error('''StartTime'' and ''StartLoc'' are mutually exlcusive')
end
if isnan(startTime) && isnan(startLoc)
    startLoc = 0;
elseif ~isnan(startTime)
    startLoc = startTime/(time(2)-time(1));
end

% Which data set are seeing(index aa) and which is the "other" data set
% (index bb)
aa = channel;

% Initialize figure
figure()
tiledlayout(sz(2),1)

% Labelling data on this axes
ax = nexttile();
hold on
plot(time, data(:,aa), 'Color', params(aa).color);
plot(time(peaks(aa).locs), data(peaks(aa).locs, aa), 'o', 'MarkerSize', 6, 'Color', peakColor);
plot(time, thresh_curve(:,aa), 'Color', threshColor)
title(params(aa).name, 'Interpreter', 'none')
xlabel('Time (s)')
ylabel(['Signal (' params(aa).units ')']);
hold off

for bb = 1:sz(2)
    % Skip the channel we are currently viewing
    if aa == bb; continue; end
    
    % This axes are available for comparision
    ax = [ax nexttile()];
    hold on
    plot(time, data(:,bb), 'Color', params(bb).color);
    plot(time(peaks(bb).locs), data(peaks(bb).locs, bb), 'o', 'MarkerSize', 6, 'Color', peakColor);
    plot(time, thresh_curve(:,bb), 'Color', threshColor)
    title(params(bb).name, 'Interpreter', 'none')
    xlabel('Time (s)')
    ylabel(['Signal (' params(bb).units ')']);
    hold off
end

start = find(peaks(aa).locs >= startLoc, 1, 'first');
numPeaks = peaks(aa).count - start + 1;

disp('Press + (=) to zoom in, - to zoom out, and r to reset zoom')
disp('Press x to cancel')
fprintf('Viewing %g peaks\n', numPeaks)
% Was the input x?
xBool = false;
marker = [];
for ii = start:peaks(aa).count
    
    pk = peaks(aa).pks(ii);
    loc = peaks(aa).locs(ii);
    
    % Mark the peak we are currently viewing
    hold(ax(1), 'on')
    delete(marker)
    marker = plot(ax(1), time(loc), pk,'x', 'MarkerSize', 8, 'Color', markerColor);
    hold(ax(1), 'off')
    
    % Determine x/y limits of the view for the current peak
    xlm = time(loc) + [-5 5];
    ylm = 1.25 * [-thresh_curve(loc,aa) pk];
    
    % Set view of current peak
    for xx = ax
        xlim(xx, xlm)
        ylim(xx, ylm)
    end
    
    % Take user input
    fprintf('Peak %g/%g\n', ii-start+1, numPeaks)
    while true
        in = input('Press enter to continue ', 's');
        switch in
            % Categorize peak based on input
            case ''
                break
                % Zoom out
            case '-'
                % Determine y limits of view for current peak
                newYlm = ylim(ax(1));
                newYlm(1) = 1.5*newYlm(1);
                newYlm(2) = 2*newYlm(2);
                
                % Set view of current peak
                for xx = ax
                    ylim(xx, newYlm)
                end
                
                % Zoom in
            case '='
                % Determine x/y limits of the view for the current peak
                newXlm = time(loc) + [-2 2];
                newYlm = [-thresh_curve(loc,aa) 1.05*pk];
                
                % Set view of current peak
                for xx = ax
                    xlim(xx, newXlm)
                    ylim(xx, newYlm)
                end
                
                % Reset view
            case 'r'
                for xx = ax
                    xlim(xx, xlm)
                    ylim(xx, ylm)
                end
                % Cancel labelling session
            case 'x'
                xBool = true;
                break
                
                % Take input until the user inputs a valid input
            otherwise
                disp('Not a valid input')
        end
    end
    if xBool; break; end
end
disp('-------------------------------------------------------------------')
end