function DiFC_Blue_maFiltData_Feb2019(src, event, disptime, window, params, channels)
%
% Display plot of moving-average filtered data in real time.
%
% V.Pera, 15 June 2016
% X. Tan, Nov, 2016; add two channels;Feb. 2019, edit for use on voltage PMTs


persistent myTime myData yyF h2

% filter
filtIdx = numel(myTime) - (window-2);

myTime = [myTime; event.TimeStamps];
myData = [myData; event.Data];


if filtIdx > 0
    
    yF1 = filter(ones(1,window)/window, 1, myData(filtIdx:end,1));
    yF2 = filter(ones(1,window)/window, 1, myData(filtIdx:end,2));
    
    yF3 = filter(ones(1,window)/window, 1, myData(filtIdx:end,3));
    yF4 = filter(ones(1,window)/window, 1, myData(filtIdx:end,4));
    % trim
    yF1 = yF1(end-src.NotifyWhenDataAvailableExceeds+1:end);
    yF2 = yF2(end-src.NotifyWhenDataAvailableExceeds+1:end);
    
    yF3 = yF3(end-src.NotifyWhenDataAvailableExceeds+1:end);
    yF4 = yF4(end-src.NotifyWhenDataAvailableExceeds+1:end);
    
else
    yF1 = filter(ones(1,window)/window, 1, myData(:,1));
    yF2 = filter(ones(1,window)/window, 1, myData(:,2));
    
    yF3 = filter(ones(1,window)/window, 1, myData(:,3));
    yF4 = filter(ones(1,window)/window, 1, myData(:,4));
end

if ~isempty(params)% if not empty
    yF1 = yF1 * params.igain1;
    yF2 = yF2 * params.igain2;
    
    yF3 = yF3 * params.igain3;
    yF4 = yF4 * params.igain4;
    myYlabel = 'Signal (mV)';
else
    myYlabel = 'Signal (V)';
end

if isempty(yyF)
    yyF = [yF1, yF2, yF3, yF4];
else
    yyF = [yyF(:,1),yyF(:,2),yyF(:,3),yyF(:,4); yF1,yF2,yF3,yF4];
end


% plot
startIdx = numel(myTime) - src.Rate*disptime+1;

if startIdx > 0 % trim data displayed according to "disptime"
    yyF = yyF(startIdx:end,:);
    myTime = myTime(startIdx:end);
    myData = myData(startIdx:end,:);
end

if isempty(h2)
    h2 = figure;
end


figure(h2);
subplot(4,1,1)
plot(myTime, yyF(:,1),'Color', 'b')
ylabel(myYlabel);

miny = min(yyF(:,1)) - 10;
maxy = max(yyF(:,1)) + 10;
axis([myTime(1) myTime(end) miny maxy])
title([channels(1).name, ': Filtered Data'])

subplot(4,1,2)
plot(myTime, yyF(:,2),'Color', 'r')
ylabel(myYlabel);
xlabel('Time (s)')
miny = min(yyF(:,2)) - 10;
maxy = max(yyF(:,2)) + 10;
axis([myTime(1) myTime(end) miny maxy])
title([channels(2).name, ': Filtered Data'])

subplot(4,1,3)
plot(myTime, yyF(:,3),'Color', 'b')
ylabel(myYlabel);
xlabel('Time (s)')
miny = min(yyF(:,3)) - 10;
maxy = max(yyF(:,3)) + 10;
axis([myTime(1) myTime(end) miny maxy])
title([channels(3).name, ': Filtered Data'])

subplot(4,1,4)
plot(myTime, yyF(:,4),'Color','r')
ylabel(myYlabel);
xlabel('Time (s)')
miny = min(yyF(:,4)) - 10;
maxy = max(yyF(:,4)) + 10;
axis([myTime(1) myTime(end) miny maxy])
title([channels(4).name, ': Filtered Data'])