function DiFC_Blue_ABsaveDataFilt_Feb2019(src, event, fname, saveInterval, params, window, disptime, channels)

% V.Pera,X.Tan 15 June 2016
% X. Tan, Nov, 2016; add two channels; Feb,2019 edit for use on voltage
% PMTS

persistent fcount scount mydata time h1;

if isempty(fcount)
    fcount = 1;
end

if isempty(scount)
    scount = 0;
end


% accumulate data
time = [time; event.TimeStamps];
mydata = [mydata; event.Data];


% filter and plot data in real time

DiFC_Blue_maFiltData_Feb2019(src, event, disptime, window, params, channels)


% save and plot data
scount = scount + 1;

if scount == saveInterval
    
    % save data
    myfile1 = [fname,'_F1','_',num2str(fcount),'.mat'];
    if isempty(params)
        data1(:,1) = mydata(:,1);
        data1(:,2) = mydata(:,2);
        save(myfile1, 'data', 'time')
        myYlabel1 = 'Signal (V)';
    else
        data1(:,1) = mydata(:,1) * params.igain1;
        data1(:,2) = mydata(:,2) * params.igain2;
        data=data1;
        save(myfile1, 'data', 'time', 'params')
        myYlabel1 = 'Signal (mV)';
    end

    myfile2 = [fname,'_F2','_',num2str(fcount),'.mat'];
    if isempty(params)
        data2(:,1) = mydata(:,3);
        data2(:,2) = mydata(:,4);
        save(myfile2, 'data', 'time')
        myYlabel2 = 'Signal (V)';
    else
        data2(:,1) = mydata(:,3) * params.igain3;
        data2(:,2) = mydata(:,4) * params.igain4;
        data=data2;
        save(myfile2, 'data', 'time', 'params')
        myYlabel2 = 'Signal (mV)';
    end
    
    % plot results
    if isempty(h1)
        h1 = figure;  
    end
    figure(h1); 
    subplot(2,1,1)
    plot(time, data1)
    ylabel(myYlabel1);
    xlabel('Time (s)')
    %legend('Channel 0', 'Channel 1')
    legend(channels(1).name, channels(2).name)
    title('Data')
    kk = strfind(fname,'\');% find the location of '\' in fname..?? several \?
    set(h1,'name',[fname(kk(end)+1:end),'_',num2str(fcount)])

    subplot(2,1,2)
    plot(time, data2)
    ylabel(myYlabel2);
    xlabel('Time (s)')
    legend(channels(3).name, channels(4).name)
    title('Data')
    kk = strfind(fname,'\');
    set(h1,'name',[fname(kk(end)+1:end),'_',num2str(fcount)])
    
    % save plot
%     myfig = [fname,'_',num2str(fcount),'.fig'];
%     hgsave(h1, myfig)

    fcount = fcount + 1; % increment counter
    scount = 0; % reset seconds counter
    
    mydata = [];
    time = [];
    
end