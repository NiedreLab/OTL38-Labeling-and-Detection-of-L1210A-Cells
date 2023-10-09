%PreProcessing for 2 Probe 2 PMT DiFC
%Josh Pace and Max Kellish 20220526 adding NIR motion processing  

warning('off','MATLAB:legend:IgnoringExtraEntries');
clear; close all;

%**********************************************************************
%CHECK THESE BEFORE RUNNING!!!!!!
rel_thresh = [5 5];
%For coincident peak detection!!!!!! 1 is used for phanntom, 0 is used for mouse
phantom = 0;

probe_distance = 3;
%**********************************************************************



%Determing if processing on a windows or mac
if ispc
    slash = '\';
else
    slash = '/';
end

noProc = 0;
plotFlag = 0;

if phantom
    dataDir = uigetdir("");
else
    dataDir = uigetdir("");
end
if dataDir == 0
    fprintf('\nProcessing canceled\n\n')
    return
end
saveDirName = sprintf([slash 'Processed_Data_and_Figures_relThresh_%g_%g'],rel_thresh(1), rel_thresh(2));
saveDir = [dataDir saveDirName];
if ~exist(saveDir,'dir')
    mkdir(saveDir)
else
    plotFlag = 0;
end

folderName = sprintf([slash 'relThresh_%g_%g_Figures'],rel_thresh(1), rel_thresh(2));
% folderName = sprintf('/relThresh_%g_%g_Figures',rel_thresh(1), rel_thresh(2));
figFolder = [saveDir folderName];
if ~exist(figFolder,'dir')
    mkdir(figFolder)
end
%proccodes library (where processing functions are stored)
%Add proccodes path for your own computer before running
%***************************************************

proccodes_library ="proccodes";

%***************************************************
addpath(genpath(proccodes_library));

%Getting file name
index = find(dataDir == slash,1,'last');
stem = dataDir(index+1:end);
if dataDir(end) ~= slash; dataDir = [dataDir slash]; end
fname = [dataDir stem];
processedDataSave = [saveDir slash stem];

procFile = sprintf('%s_proc_relThresh_%g_%g.mat',processedDataSave, rel_thresh(1), rel_thresh(2));
if exist(procFile, 'File')
    noProc = 1;
    plotFlag = 0;
    load(procFile);
end
%-----------------------------------------------------------------------------------------------------------------%
%% Running proccessing code
fprintf('Running 2 Probe 2 PMT preProc for %s...\n\n',stem)
if ~noProc
    fprintf('Pre-processing %s...\n',stem)
    plotFlag = 1;
    runName = input("Scan name for plot title: ",'s');
    titleName(1).name = strcat(runName, " P1");
    titleName(2).name = strcat(runName, " P2");
    % Load data from both fibers
    load([fname '_F1.mat'], 'time', 'data', 'params')
    data_1 = data(:,1);
    params1 = params;
    load([fname '_F2.mat'], 'data', 'params')
    %- sign for when probe 1 is lock In and probe 2 is not
    data_2 = data(:,1);
    params2 = params;
    
    % Sampling frequnecy
    fs = 1 ./ (time(2) - time(1));

    % Format data from both fibers
    data = [data_1 data_2];
    sources = length(data(1,:));
    params = [params1(1) params2(1)];
    params(1).name = [stem ' Probe 1'];
    params(2).name = [stem ' Probe 2'];
    clear data_1 data_2 params1 params2
    if params(1).units == 'mV'
        systemName = 'BG';
    else
        systemName = 'NIR';
    end

%%

    % Pre-proc background subtracts data, calculates the noise/ moving peak
    % threshhold, and identifies peak candidates. Processes all data in the
    % data array
    if phantom
        [data_bs, noise, peaks, thresh_curve,in_dat] = preProc(data, time,[], 'RelativeThresh', rel_thresh,'RemoveBunchedPeaks','False');
    else
        [data_bs, noise, peaks, thresh_curve,in_dat] = preProc(data, time,[], 'RelativeThresh', rel_thresh);
    end
    %,'Mode','setThresh', 'HardThresh', []
    % Calculating SNR based on the estimated run noise, not control noise
    
    %Getting all Peaks SNR
    for ii = 1:sources
        all_snr(ii).dB = 20*log10(peaks(ii).pks./noise(ii));
    end
    
    %Getting peaks/min
    for ii = 1:sources
        pk_per_min(ii).num = peaks(ii).count/time(end)*60;
    end
    %Getting detection time for the count rate plot
    for ii = 1:sources
        peakTime(ii).time = time(peaks(ii).locs);
    end
    if length(peakTime(1).time) >= length(peakTime(2).time)
        detection_time = sort(peakTime(1).time);
    else
        detection_time = sort(peakTime(2).time);
    end
    
    %Time of scan in seconds
    scan_length = length(time)/fs;
    %Calculating mean background level of the orginal data
    mean_background = mean(-data);
    
    %-----------------------------------------------------------------------------------------------------------------%
    % REMOVE COINCIDENT PEAKS
    %-----------------------------------------------------------------------------------------------------------------%
    %Used for in vivo scans, not phantom
    if phantom
        coinc_pk_count = 0;
        coinc_peaks = 0;
    else
    disp('Looking for any coincident peaks...');
   disp('Looking for any coincident peaks...');
    [coinc_peaks(1), coinc_peaks(2),scores] = matchCoincPeaks(peaks(1), peaks(2), time);
    [peaks(1), peaks(2), coinc_pk_count] =... 
    removeCoincPeaks(peaks(1), peaks(2), time, 'CoincidenceWindow',0.03);
    end
    %plotPeaks(data_bs, time,peaks, thresh_curve, params);
    %-----------------------------------------------------------------------------------------------------------------%
    % MATCH PEAKS IN DIRECTIONS
    %-----------------------------------------------------------------------------------------------------------------%
    disp('Matching peaks in the forward direction...')
  [fwd_peaks(1), fwd_peaks(2), fwd_speed, fwd_score] =...
        matchDirectionalPeaks(peaks(1), peaks(2), time, probe_distance);
    %plotPeaks(data_bs, time, fwd_peaks, thresh_curve, params,'Direction', 'fwd');
     

   
    disp('Matching peaks in the reverse direction...')
    [rev_peaks(1), rev_peaks(2), rev_speed, rev_score] =...
        matchDirectionalPeaks(peaks(1), peaks(2), time, probe_distance,'Direction','rev');
    %plotPeaks(data_bs, time, rev_peaks, thresh_curve, params,'Direction', 'rev');

%     %Looking at potential double matches in Probe 1
    [fwd_peaks, rev_peaks, fwd_score, rev_score, ...
     fwd_speed,rev_speed, doubleMatchCountsP1] = ...
     determineBestDirectionalMatch(fwd_peaks, rev_peaks, fwd_score, rev_score,fwd_speed,rev_speed,1);
     fprintf("Found %g double matches when looking at Probe 1\n",doubleMatchCountsP1);
     %Looking at potential double matches in Probe 2
     [fwd_peaks, rev_peaks, fwd_score, rev_score, ...
     fwd_speed,rev_speed, doubleMatchCountsP2] = ...
     determineBestDirectionalMatch(fwd_peaks, rev_peaks, fwd_score, rev_score,fwd_speed,rev_speed,2);
     fprintf("Found %g double matches when looking at Probe 2\n",doubleMatchCountsP2);
         %Getting fwd matched peak SNR
    for ii = 1:sources
        fwd_match_snr(ii).dB = 20*log10(fwd_peaks(ii).pks./noise(ii));
    end

     %Getting reverse matched peak SNR
    for ii = 1:sources
        rev_match_snr(ii).dB = 20*log10(rev_peaks(ii).pks./noise(ii));
    end
    
    %Getting a struct with only the unmatched peaks
    unmatched_peaks = peaks;
    for ii = 1:2
        % build unmatched peak struct
        unmatched_peaks_arr = true(peaks(ii).count, 1);
        
        [~, jj, ~] = intersect(peaks(ii).locs, fwd_peaks(ii).locs);
        unmatched_peaks_arr(jj) = false;
        
        [~, jj, ~] = intersect(peaks(ii).locs, rev_peaks(ii).locs);
        unmatched_peaks_arr(jj) = false;
        
        unmatched_peaks(ii).pks = unmatched_peaks(ii).pks(unmatched_peaks_arr);
        unmatched_peaks(ii).locs = unmatched_peaks(ii).locs(unmatched_peaks_arr);
        unmatched_peaks(ii).widths = unmatched_peaks(ii).widths(unmatched_peaks_arr);
        unmatched_peaks(ii).proms = unmatched_peaks(ii).proms(unmatched_peaks_arr);
        
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
     
   
    
        % Save processed data, including parameters for processing

    pname = sprintf('%s_proc_relThresh_%g_%g.mat', processedDataSave, rel_thresh(1), rel_thresh(2));
    save(pname, 'time','data' ,'data_bs', 'params', 'noise','all_snr','fwd_match_snr','rev_match_snr',...
             'peaks','unmatched_peaks','fwd_peaks','rev_peaks','coinc_peaks','fwd_score', 'fwd_speed','rev_score','rev_speed',...
             'mean_background', 'thresh_curve','pk_per_min','coinc_pk_count','sources',...
             'detection_time','scan_length', 'in_dat');
    disp('Saved processed data...')
    
    % Printing information after preprocessing 
    for ii = 1:sources
        fprintf('Detected %4.2g peaks per min in Probe %g\n', pk_per_min(ii).num,ii)
    end
    for ii = 1:sources
        fprintf('Peaks candidates found for Probe %g: %g\n', ii,peaks(ii).count);
    end
    for ii = 1:sources
        fprintf('Calculated noise for Probe %g: %.3f %s\n',ii,noise(ii),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Estimated All Peak Probe %g SNR: %.3f dB\n',ii, mean(all_snr(ii).dB));
    end
    for ii = 1:sources
        fprintf('Estimated Fwd Matched Probe %g SNR: %.3f dB\n',ii, mean(fwd_match_snr(ii).dB));
    end
    for ii = 1:sources
        fprintf('Estimated Rev Matched Probe %g SNR: %.3f dB\n',ii, mean(rev_match_snr(ii).dB));
    end
    for ii = 1:sources
        fprintf('Average relative threshold for Probe %g: %.3f %s\n', ii,mean(thresh_curve(ii)),params(ii).units);
    end
    for ii = 1:sources
       fprintf('Average background level for Probe %g: %.3f %s\n',ii,mean_background(ii),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Average peak amplitude for Probe %g: %.3f %s\n',ii, mean(peaks(ii).pks),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Peak amplitude standard deviation for Probe %g: %.3f %s\n',ii, std(peaks(ii).pks),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Largest peak amplitude for Probe %g: %.3f %s\n', ii,max(peaks(ii).pks),params(ii).units);
        fprintf('Smallest peak amplitude for Probe %g: %.3f %s\n', ii,min(peaks(ii).pks),params(ii).units);
    end
    
    try
        fprintf('Removed a total of %g coincident peaks\n',coinc_pk_count);
        fprintf('Arterial Direction Matched Peaks: %g\n', fwd_peaks(1).count);
        fprintf('Venous Direction Matched Peaks: %g\n', rev_peaks(1).count);
    catch
    end
else
    % Printing the infromation if the data has already been processed
    for ii = 1:sources
        fprintf('Detected %4.2g peaks per min in Probe %g\n', pk_per_min(ii).num,ii)
    end
    for ii = 1:sources
        fprintf('Peaks candidates found for Probe %g: %g\n', ii,peaks(ii).count);
    end
    for ii = 1:sources
        fprintf('Calculated noise for Probe %g: %.3f %s\n',ii,noise(ii),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Estimated All Peak Probe %g SNR: %.3f dB\n',ii, mean(all_snr(ii).dB));
    end
    for ii = 1:sources
        fprintf('Estimated Fwd Matched Probe %g SNR: %.3f dB\n',ii, mean(fwd_match_snr(ii).dB));
    end
    for ii = 1:sources
        fprintf('Estimated Rev Matched Probe %g SNR: %.3f dB\n',ii, mean(rev_match_snr(ii).dB));
    end
    for ii = 1:sources
        fprintf('Average relative threshold for Probe %g: %.3f %s\n', ii,mean(thresh_curve(ii)),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Average background level for Probe %g: %.3f %s\n',ii,mean_background(ii),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Average peak amplitude for Probe %g: %.3f %s\n',ii, mean(peaks(ii).pks),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Peak amplitude standard deviation for Probe %g: %.3f %s\n',ii, std(peaks(ii).pks),params(ii).units);
    end
    for ii = 1:sources
        fprintf('Largest peak amplitude for Probe %g: %.3f %s\n', ii,max(peaks(ii).pks),params(ii).units);
        fprintf('Smallest peak amplitude for Probe %g: %.3f %s\n', ii,min(peaks(ii).pks),params(ii).units);
    end
    
    try
        fprintf('Removed a total of %g coincident peaks\n',coinc_pk_count);
        fprintf('Arterial Direction Matched Peaks: %g\n', fwd_peaks(1).count);
        fprintf('Venous Direction Matched Peaks: %g\n', rev_peaks(1).count);
    catch
        fprintf('Removed a total of %g coincident peaks\n',coinc_pk_count);
        fprintf('Arterial Direction Matched Peaks: %g\n', fwd_peaks(1).count);
    end
end

if plotFlag
    disp("Plotting and saving plots...");
    %PLAYING WITH PLOTS
    lgnd = {params.name};
    load([fname '_F1.mat'], 'data', 'time')
    originalTime = time;
    data_1 = data(:,1);
    load([fname '_F2.mat'], 'data')
    data_2 = data(:,1);
    data = [data_1  data_2];
    sources = length(data(1,:));
    clear data_1 data_2
    
    if exist('filteredData', 'var')
%         %Cutting out first 10 seconds because of the filter behavior around 0
    timeInSeconds = 5;
    timeShift = (timeInSeconds*2000)+1;%How much time to cut off from the beginning because of the filtering
    filterData_1 = data_bs(timeShift:end,1);
    filterData_2 = data_bs(timeShift:end,2);
    data_bs = [filterData_1 filterData_2];
    filterThreshCurve_1 = thresh_curve(timeShift:end,1);
    filterThreshCurve_2 = thresh_curve(timeShift:end,2);
    thresh_curve = [filterThreshCurve_1 filterThreshCurve_2];
    filterTime = time(timeShift:end);
    time = filterTime-timeInSeconds;

        truePeakCandidates_1 = peaks(1).locs > timeShift;
        truePeakCandidates_2 = peaks(2).locs > timeShift;
        peaks(1) = filterPeaks(peaks(1), truePeakCandidates_1);
        peaks(2) = filterPeaks(peaks(2), truePeakCandidates_2);
      
        
        trueFwdPeakCandidates_1 = fwd_peaks(1).locs > timeShift;
        trueFwdPeakCandidates_2 = fwd_peaks(2).locs > timeShift;
        fwd_peaks(1) = filterPeaks(fwd_peaks(1), trueFwdPeakCandidates_1);
        fwd_peaks(2) = filterPeaks(fwd_peaks(2), trueFwdPeakCandidates_2);
        
        
        
        trueRevPeakCandidates_1 = rev_peaks(1).locs > 2001;
        trueRevPeakCandidates_2 = rev_peaks(2).locs > 2001;
        
        rev_peaks(1) = filterPeaks(rev_peaks(1), trueRevPeakCandidates_1);
        rev_peaks(2) = filterPeaks(rev_peaks(2), trueRevPeakCandidates_2);
        
        
        trueCoincidentPeakCandidates_1 = coinc_peaks(1).locs > timeShift;
        trueCoincidentPeakCandidates_2 = coinc_peaks(2).locs > timeShift;
        coinc_peaks(1) = filterPeaks(coinc_peaks(1), trueCoincidentPeakCandidates_1);
        coinc_peaks(2) = filterPeaks(coinc_peaks(2), trueCoincidentPeakCandidates_2);
        
    for ii = 1:sources    
%         for i = 1:length(peaks(1).locs)
%             if peaks(ii).locs(i) < 20001
%                 peaks(ii).locs(i) = [];
%             end
%         end
        peaks(ii).locs = peaks(ii).locs-timeShift;
        
        fwd_peaks(ii).locs = fwd_peaks(ii).locs-timeShift;
        
        rev_peaks(ii).locs = rev_peaks(ii).locs-timeShift;
        
        coinc_peaks(ii).locs = coinc_peaks(ii).locs-timeShift;
        
        
        
        for i = 1:length(peaks(ii).locs)
            if peaks(ii).locs(i) < 0
                peaks(ii).locs(i) = -1*peaks(ii).locs(i);
            end
        end
        
        for i = 1:length(fwd_peaks(ii).locs)
            if fwd_peaks(ii).locs(i) < 0
                fwd_peaks(ii).locs(i) = -1*fwd_peaks(ii).locs(i);
            end
        end
        
        for i = 1:length(rev_peaks(ii).locs)
            if rev_peaks(ii).locs(i) < 0
                rev_peaks(ii).locs(i) = -1*rev_peaks(ii).locs(i);
            end
        end
        
        for i = 1:length(coinc_peaks(ii).locs)
            if coinc_peaks(ii).locs(i) < 0
                coinc_peaks(ii).locs(i) = -1*coinc_peaks(ii).locs(i);
            end
        end
        
%         for i = 1:length(fwd_peaks(1).locs)
%                     if fwd_peaks(ii).locs(i) < 20001
%                         fwd_peaks(ii).locs(i) = [];
%                     end
%         end
%          for i = 1:length(rev_peaks(1).locs)
%                     if rev_peaks(ii).locs(i) < 20001
%                         rev_peaks(ii).locs(i) = [];
%                     end
%          end
% 
%          for i = 1:length(coinc_peaks(1).locs)
%                 if coinc_peaks(ii).locs(i) < 20001
%                     coinc_peaks(ii).locs(i) = [];
%                 end
%          end
    end
    end

    %Getting y axis limts so both probes have the same y axis
    dataMax = max(max(data_bs))
    dataMin = min(min(data_bs))
    
    
%     yAxesMaxTmp = round(dataMax,-2);
%     
%     if yAxesMaxTmp == 0
%         yAxesMax = 50;
%     else
%         if yAxesMaxTmp < dataMax
%             yAxesMax = 100 + yAxesMaxTmp;
%         else
%             yAxesMax = yAxesMaxTmp;
%         end
%     end
%     if le(dataMin, -50)
%         yAxesMin = -100;
%     else
%         yAxesMin = -50;
%     end
%     
%     yAxesMax  = 200;
%     yAxesMin = -200;

    yAxesMax = input("Please set the Y axis max value: ");
    yAxesMin = input("Please set the Y axis min value: ");
    
    if yAxesMin > 0
        yAxesMin = -yAxesMin;
    end
    % Plot original data
    figure(1);
    tiledlayout(2,1)
    for ii = 1:sources
        nexttile()
        plot(originalTime, data(:,ii), 'Color',params(ii).color);
        hold on;
        title([titleName(ii).name,' Original Data'], 'Interpreter', 'none','FontSize',18')
        legend('DiFC Data')
        axis = gca;
        axis.FontSize = 18;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        hold off
    end
    fig = gcf;
    tl= fig.Children;
    xlabel(tl,'Time (s)','FontSize',18,'FontWeight','bold');
    ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',18,'FontWeight','bold')
    saveas(figure(1),[figFolder slash 'Original Data.fig']); 
    % Plot smoothed, background subtracted data
    figure(2);
    tiledlayout(2,1)
    for ii = 1:sources 
        nexttile()
        plot(time, data_bs(:,ii), 'Color',params(ii).color);
        hold on;
        title([titleName(ii).name,' Smoothed, Bkgnd Subtracted Data'], 'Interpreter', 'none','FontSize',18')
        legend('DiFC Data')
        axis = gca;
        axis.FontSize = 18;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        hold off;
        ylim([yAxesMin yAxesMax])
    end
    fig = gcf;
    tl= fig.Children;
    xlabel(tl,'Time (s)','FontSize',18,'FontWeight','bold');
    ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',18,'FontWeight','bold')
%     saveas(figure(2),[figFolder slash 'Smoothed and Bkgrnd Subtracted Data.fig']); 
    % Plot Peak Candidates
    figure(3)
    tiledlayout(2,1)
    for ii = 1:sources
        nexttile()
        plot(time, data_bs(:,ii), 'Color',params(ii).color)
        hold on
        plot(time(peaks(ii).locs), data_bs(peaks(ii).locs,ii), 'mo', 'MarkerSize', 6)
        plot(time, thresh_curve(:,ii), '-k')
        title([titleName(ii).name, ' Peak Candidates'], 'Interpreter', 'none','FontSize',18')
        legend('DiFC Data', sprintf('Peaks: %g',peaks(ii).count), 'Dynamic Threshold')
        axis = gca;
        axis.FontSize = 18;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        hold off
        ylim([yAxesMin yAxesMax])
    end
    fig = gcf;
    tl= fig.Children;
    xlabel(tl,'Time (s)','FontSize',18,'FontWeight','bold');
    ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',18,'FontWeight','bold')
%     saveas(figure(3),[figFolder slash 'Peak Candidates.fig']); 
    % Plot Foward Matches
    figure(4)
    tiledlayout(2,1)
    for ii = 1:sources
        nexttile()
        plot(time, data_bs(:,ii), 'Color',params(ii).color)
        hold on
        plot(time(fwd_peaks(ii).locs), data_bs(fwd_peaks(ii).locs,ii), 'm>', 'MarkerSize', 6)
        plot(time, thresh_curve(:,ii), '-k')
        title([titleName(ii).name, ' Foward Direction Matched Peaks'], 'Interpreter', 'none','FontSize',18')
        legend('DiFC Data',sprintf( 'Foward Matched Peaks: %g',fwd_peaks(ii).count), 'Dynamic Threshold')
        axis = gca;
        axis.FontSize = 18;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        hold off
        ylim([yAxesMin yAxesMax])
    end
    fig = gcf;
    tl= fig.Children;
    xlabel(tl,'Time (s)','FontSize',17,'FontWeight','bold');
    ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',18,'FontWeight','bold')
%     saveas(figure(4),[figFolder slash 'Foward Matches.fig']); 
    % Plot Reverse Matches
        figure(5)
        tiledlayout(2,1)
        for ii = 1:sources
            nexttile()
            plot(time, data_bs(:,ii), 'Color',params(ii).color)
            hold on
            plot(time(rev_peaks(ii).locs), data_bs(rev_peaks(ii).locs,ii), 'm<', 'MarkerSize', 6)
            plot(time, thresh_curve(:,ii), '-k')
            title([titleName(ii).name, ' Reverse Direction Matched Peaks'], 'Interpreter', 'none','FontSize',17')
            legend('DiFC Data',sprintf( 'Reverse Matched Peaks: %g',rev_peaks(ii).count), 'Dynamic Threshold')
            axis = gca;
            axis.FontSize = 18;
            hold off
            axis.LineWidth = 2;
            axis.FontWeight = 'bold';
            axis.FontName = 'Arial';
            ylim([yAxesMin yAxesMax])
        end
        fig = gcf;
        tl= fig.Children;
        xlabel(tl,'Time (s)','FontSize',18,'FontWeight','bold');
        ylabel(tl, [systemName '-DiFC Signal (' params(1).units ')'], 'FontSize',18,'FontWeight','bold')
%         saveas(figure(5),[figFolder slash 'Reverse Matches.fig']); 
     %     Plot Histogram of Peak Candidate Amplitudes
        try
            %binLimit =  max([max(peaks(1).pks) max(peaks(2).pks)]) + 10;
            all_peak_candidate_SNR = [peaks(1).pks; peaks(2).pks];
            [N,peak_cadidate_edges] = histcounts(all_peak_candidate_SNR);
            figure(6)
            tiledlayout(2,1)
            for ii = 1:sources
                nexttile()
                %binWidth = (max(peaks(ii).pks) - min(peaks(ii).pks))/ceil(sqrt(length(peaks(ii).pks)));
                %histogram(peaks(ii).pks, 'BinWidth', binWidth, 'BinLimit', [0 binLimit])
                histogram(peaks(ii).pks,'BinEdges',peak_cadidate_edges)
                title([titleName(ii).name, ' Peak Candidate Amplitudes'], 'Interpreter', 'none','FontSize',18')
%                 xlabel(sprintf('Peak Candidate Amplitude (%s)',params(ii).units));
%                 ylabel('Count','FontSize',15);
                axis = gca;
                axis.FontSize = 18;
                axis.LineWidth = 2;
                axis.FontWeight = 'bold';
                axis.FontName = 'Arial';
                hold off
            end
            fig = gcf;
            tl= fig.Children;
            xlabel(tl,sprintf('Peak Candidate Amplitude (%s)',params(1).units),'FontSize',18,'FontWeight','bold');
            ylabel(tl,'Count','FontSize',18,'FontWeight','bold');
            saveas(figure(6),[figFolder slash 'Peak Candidate Amplitude Histogram.fig']); 
        catch
           fprintf("No Peak Candidate Histogram displayed for Probe 2\n");
           saveas(figure(6),[figFolder slash 'Peak Candidate Amplitude Histogram.fig']); 
        end
        
        %     Plot Histogram of Fwd Peak Candidate SNRs
        try
            %binLimit =  max([max(fwd_match_snr(1).dB) max(fwd_match_snr(2).dB)]) + 10;
            all_fwd_snr = [fwd_match_snr(1).dB; fwd_match_snr(2).dB];
            [N, fwd_match_snr_edges] = histcounts(all_fwd_snr);
            figure(7)
            tiledlayout(2,1)
            for ii = 1:sources
                nexttile()
                %binWidth = (max(fwd_match_snr(ii).dB) - min(fwd_match_snr(ii).dB))/ceil(sqrt(length(fwd_match_snr(ii).dB)));
                %histogram(fwd_match_snr(ii).dB, 'BinWidth', binWidth, 'BinLimit', [0 binLimit])
                histogram(fwd_match_snr(ii).dB, 'BinEdges',fwd_match_snr_edges)
                title([titleName(ii).name, ' Foward Matched Peak Candidate SNR'], 'Interpreter', 'none','FontSize',18')
                axis = gca;
                axis.FontSize = 18;
                axis.LineWidth = 2;
                axis.FontWeight = 'bold';
                axis.FontName = 'Arial';
                hold off
            end
            fig = gcf;
            tl= fig.Children;
            xlabel(tl,'SNR (dB)','FontSize',18,'FontWeight','bold');
            ylabel(tl,'Count','FontSize',18,'FontWeight','bold');
            saveas(figure(7),[figFolder slash 'Fwd Matched Peak Candidate SNR Histogram.fig']); 
        catch
           fprintf("No Fwd Matched SNR Histogram displayed\n");
           saveas(figure(7),[figFolder slash 'Fwd Matched Peak Candidate SNR Histogram.fig']); 
        end
        %     Plot Histogram of Rev Peak Candidate SNRs
        
        try
            %binLimit =  max([max(rev_match_snr(1).dB) max(rev_match_snr(2).dB)]) + 10;
            all_rev_snr = [rev_match_snr(1).dB; rev_match_snr(2).dB];
            [N, rev_match_snr_edges] = histcounts(all_rev_snr);
            figure(8)
            tiledlayout(2,1)
            for ii = 1:sources
                nexttile()
                %binWidth = (max(rev_match_snr(ii).dB) - min(rev_match_snr(ii).dB))/ceil(sqrt(length(rev_match_snr(ii).dB)));
                %histogram(rev_match_snr(ii).dB, 'BinWidth', binWidth, 'BinLimit', [0 binLimit])
                histogram(rev_match_snr(ii).dB, 'BinEdges', rev_match_snr_edges)
                title([titleName(ii).name, ' Reverse Matched Peak Candidate SNR'], 'Interpreter', 'none','FontSize',18')
                axis = gca;
                axis.FontSize = 18;
                axis.LineWidth = 2;
                axis.FontWeight = 'bold';
                axis.FontName = 'Arial';
                hold off
            end
            fig = gcf;
            tl= fig.Children;
            xlabel(tl,'SNR (dB)','FontSize',18,'FontWeight','bold');
            ylabel(tl,'Count','FontSize',18,'FontWeight','bold');
            saveas(figure(8),[figFolder slash 'Rev Matched Peak Candidate SNR Histogram.fig']); 
        catch
           fprintf("No Rev Matched SNR Histogram displayed\n");
           saveas(figure(8),[figFolder slash 'Rev Peak Candidate SNR Histogram.fig']); 
        end
        if sources >=2 
            params(1).name = titleName(1).name;
            params(2).name = titleName(2).name;
            pc = 0;
            plotAllPeaks3(data_bs, time, peaks, fwd_peaks, rev_peaks, coinc_peaks,thresh_curve, params,figFolder,phantom,yAxesMax, yAxesMin)
        end
        
        
        
        
        %Getting and plotting detection count rate
        
        
        % Set the length of the window/interval in seconds (tail artery is ~50 uL
        % per minute, 2000uL total blood volume)
        % 120s = 5% of blood vol
        interval_length = 120; %seconds
        
        time0 = linspace(0,scan_length, scan_length*fs)';

        % Scan average number of detections per minute
        avg_events = length(detection_time)./ (scan_length/60);
        x_axis_span = [0.01 scan_length/60 + 0.01];

        % Calculate sliding window
        CTCs_per_interval = countDetectionsPerInterval(detection_time, interval_length, scan_length, fs);
        time = time0(interval_length*fs/2:end-interval_length*fs/2)/60;

        % Scan average number of detections per interval
        avg = avg_events * interval_length / 60;
        
        % Plotting! 
        figure(10);
        tiledlayout(1,1)
        nexttile;
        plot(x_axis_span, [avg avg], 'k-', 'LineWidth', 1.5); 
        hold on;
        plot(time,CTCs_per_interval','b-','LineWidth', 1.5);
        title([runName, ' Detection Count Rate'], 'Interpreter', 'none','FontSize',17')
        xlabel('Time (min)','FontSize',17,'FontWeight','bold');
        ylabel('Count Rate','FontSize',17,'FontWeight','bold');
        legend('Scan Average')
        axis = gca;
        axis.FontSize = 18;
        axis.LineWidth = 2;
        axis.FontWeight = 'bold';
        axis.FontName = 'Arial';
        hold off;
        saveas(figure(10),[figFolder  slash 'Scan Count Rate.fig']); 
end

fprintf("\n\nProcessing and displaying results complete for: %s\n\n",stem);

% seePeaks(data_bs,time,peaks,thresh_curve,params);