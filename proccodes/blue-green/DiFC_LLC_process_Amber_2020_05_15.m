function [out_dat] = DiFC_LLC_process_Amber_2020_05_15(data_path, data_file_name, save_file, save_to )
% Edit of Jessica's Code DiFC_LLC_batch_process_NRT downloaded on
% 2020_04_21

%DiFC_LLC_process_one -  Data processing routine for DiFC data, specificially
%the 'voltage output' system with LLC tumor metatasis model
% M.Niedre May 2019; Revised Nov 2019
% code written to analyze DiFC data in LLC lung metastasis model.
% the code expects two files corresponding to measurements from fiber 1 and 2, respectively (_F1, _F2).

%%   User-Defined Parameters
data_file = strcat(data_path, data_file_name);
plot_flag = 1;
remove_10min =  0;                                      % flag to remove first 10 minutues of data (=1) which can be noisy

%  Parameters for data pre-processing
background_sliding_window = 2.5;                % length of median filter used for background subtraction in seconds
smoothing_siding_window = 0.005;                % length of median filter used for data filtering in seconds

%   Parameters for  peak detection routine
rel_thresh = 4.8;                                     % multiplicative factor for the adaptive threshold;  threshold is rel_thresh*estimated standard deviation
prominence_factor = 1;                            % minimum peak prominence (versus adjacent peaks) is the threshold/prominence_factor
coinc_window = 0.03;                            % time window to search for coincident peaks in s
%  Parameters for ouput
est_cell_amplitude = 15;                    % Estimated amplitude of a single LLC cell in mV; used for count rate estimation
super_peak_factor = 3;

%% Data Procesisng
%   Reads in data from the two fibers, performs pre-processing, and peak  detection
fprintf('\n\n\n\n\n\n\nDiFC Data Processing Routine for LLC Mice\n')
%loads the user-supplied raw data (data_file) file for fiber 1
dat_f   = strcat(data_file, '_F1.mat'); load(dat_f);
fprintf(strcat('\nLoading raw data file for fiber 1:\t', data_file, '\n'))

% first check that 'time' and 'data' variables exist, if they don't then notify user and break; otherwise proceed
if exist('time') ~= 1 || exist('data') ~= 1
    fprintf('Error: time and data variables not found - check raw data file\n')
    return
end
% extract basic time / sampling parameters
dt      = time(2) - time(1);  % time increment
fs      = 1 / dt;  % sampling rate used in acquisition
bsw     = background_sliding_window * fs; bsw = int16(bsw); % background sliding window
ssw     = smoothing_siding_window * fs;    % smoothing sliding window

%check for noisy/bad PMT measurement (occasionally happens)
fprintf('\tChecking for bad/noisy PMTs on Fiber 1\n')
[data] = check_noisy_pmt(data, bsw, 1);

% pre-process the sum of channel 1 and 2 for fiber 1.  Finds peak
% candidates based on a calculated adaptive threshold.
fprintf('\tPre-processing data for Fiber 1 and looking for peaks\n')
data_sum = data(:,1)+data(:,2); if sign(data_sum) ==-1; data_sum =-data_sum; end
if remove_10min == 1
    fprintf('\tWarning: First 10 minutes removed from analysis.  If not intended clear flag "remove_10min" \n')
    data_sum = data_sum(600/dt+1:end);
    time = time(1:end-600/dt);
end
data_sum = smooth(data_sum, ssw);
[ch_pp(:,1), ch_sd(1), ch_bg(:,1), ch1_info]  = pre_proc_ch(data_sum, fs, bsw, time, rel_thresh, prominence_factor);
clear data;


%load and processes the raw data for fiber 2 and find peak candidates (repeat of the above steps for the second data file)
dat_f   = strcat(data_file, '_F2.mat'); load(dat_f);
fprintf(strcat('\nLoading raw data file for fiber 2:\t', data_file, '\n'))
if exist('time') ~= 1 || exist('data') ~= 1
    fprintf('Error: time and data variables not found - check raw data file\n')
    return
end

%check for noisy/bad PMT measurement (occasionally happens)
fprintf('\tChecking for bad/noisy PMTs on Fiber 2\n')
[data] = check_noisy_pmt(data, bsw ,2);
fprintf('\tPre-processing data for Fiber 2 and looking for peaks\n')
data_sum = data(:,1)+data(:,2); if sign(data_sum) ==-1; data_sum =-data_sum; end
if remove_10min == 1
    data_sum = data_sum(600/dt+1:end);
    time = time(1:end-600/dt);
end
data_sum = smooth(data_sum, ssw);
[ch_pp(:,2), ch_sd(2) ch_bg(:,2), ch2_info]  = pre_proc_ch(data_sum, fs, bsw, time, rel_thresh, prominence_factor);
clear data;

%%
fprintf('\n\n')
scan_length_min = time(end)/60;
fprintf(strcat('Scan Length in Minutes:\t\t\t\t', num2str(scan_length_min),'\n'));

fprintf(strcat('Channel 1 Mean estimated standard deviation:\t', num2str(ch_sd(1)),'\n'));
fprintf(strcat('Channel 1 Peak candidates identified:\t\t', num2str(length(ch1_info.locs)),'\n'));

fprintf(strcat('Channel 2 Mean Estimated Standard Deviation:\t', num2str(ch_sd(2)),'\n'));
fprintf(strcat('Channel 2 Peak candidates identified:\t\t', num2str(length(ch2_info.locs)),'\n'));

% generate a plot of the background-subtracted signal, peak candidates,
% and the calculated adaptive threshold.
if plot_flag == 1
    figure;
    %subplot(2,2,1); 
    %plot(time, ch_pp(:,1)); hold on;
    ch1_info.locs = int32(ch1_info.locs);
    ch2_info.locs = int32(ch2_info.locs);
    plot(time(ch1_info.locs), ch_pp(ch1_info.locs,1), 'mo', 'MarkerSize', 4);
    plot(time, ch1_info.thresh_curve, 'g-'); xlabel('Time (s)'); ylabel('DiFC Signal (mV)');
  %subplot(2,2,3);
   plot(time, ch_pp(:,2)); hold on;
    
   plot(time(ch2_info.locs), ch_pp(ch2_info.locs,2), 'mo', 'MarkerSize', 4 ); plot(time, ch2_info.thresh_curve, 'g-'); xlabel('Time (s)'); ylabel('DiFC Signal (mV)');
end

% Search for coincident peaks indicative of movement artifacts
fprintf('\n\nSearching for coincident peaks (motion artifacts)\n')
[ch1_info, ch2_info, coinc_pk_count, multi_pk_count] = find_coinc_peaks(ch1_info, ch2_info, fs, coinc_window);
if(coinc_pk_count ~=0)
    fprintf(strcat('   Warning: ', num2str(coinc_pk_count), ' coincident peaks identified and removed from the peak count\n'))
end
if(multi_pk_count ~=0)
    fprintf(strcat('   Warning:   ', num2str(multi_pk_count), ' incidences of grouped peaks identified and removed from the peak count\n'))
end

if plot_flag == 1
    figure; subplot(2,2,1); hold on; plot(time(ch1_info.locs), ch_pp(ch1_info.locs,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); legend('Data', 'Peak Candidates', 'Adaptive Threshold', 'Corrected Peaks');
    figure; subplot(2,2,3); hold on; plot(time(ch2_info.locs), ch_pp(ch2_info.locs,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); legend('Data', 'Peak Candidates', 'Adaptive Threshold', 'Corrected Peaks' )
end

cor_pks1 = length(ch1_info.pks);
cor_pks2 = length(ch2_info.pks);
est_cells1 = round(ch1_info.pks./est_cell_amplitude); tot_est_cells1 = sum(est_cells1);
est_cells2 = round(ch2_info.pks./est_cell_amplitude); tot_est_cells2 = sum(est_cells2);
super_peaks1 = find(est_cells1>=super_peak_factor); num_super_peaks1 = length(super_peaks1); %if isempty(super_peaks1); super_peaks1 = 0; end
super_peaks2 = find(est_cells2>=super_peak_factor); num_super_peaks2 = length(super_peaks2); %if isempty(super_peaks2); super_peaks2 = 0; end

if plot_flag == 1
    figure;
    subplot(2,2,2);
    if cor_pks1 >=1
        for j = 1:cor_pks1; plot(time(ch1_info.locs(j))/60, .25, 'bo', 'MarkerSize', min(50,(6+round(est_cells1(j)))), 'LineWidth', 3); hold on;  end
    end
    ylim([0 .5]); xlim([0 max(time)/60]); xlabel('Time (minutes)'); ylabel('Detections - Ch1');
    subplot(2,2,4);
    if cor_pks2 >=1
        for j = 1:cor_pks2;  plot(time(ch2_info.locs(j))/60, 0.25, 'ro', 'MarkerSize', min(50,(6+round(est_cells2(j)))), 'LineWidth', 3); hold on;  end
    end
    ylim([0 .5]); xlim([0 max(time)/60]); xlabel('Time (minutes)'); ylabel('Detections - Ch2');
end
fprintf(strcat('\n\nChannel 1 corrected peak count:\t\t\t', num2str(cor_pks1),'\n'));
fprintf(strcat('Channel 2 corrected peak count:\t\t\t', num2str(cor_pks2),'\n'));
fprintf(strcat('\nChannel 1 corrected peak count rate (/min):\t', num2str(cor_pks1./scan_length_min),'\n'));
fprintf(strcat('Channel 2 corrected peak count rate (/min):\t', num2str(cor_pks2./scan_length_min),'\n'));
fprintf(strcat('\nChannel 1 estimated cell detection rate (/min):\t', num2str(tot_est_cells1./scan_length_min),'\n'));
fprintf(strcat('Channel 2 estimated cell detection rate (/min):\t', num2str(tot_est_cells2./scan_length_min),'\n'));
fprintf(strcat('\nTotal scan peak count rate (/min):\t\t', num2str((cor_pks1+cor_pks2)./scan_length_min),'\n'));
fprintf(strcat('Total  estimated cell detection rate (/min):\t', num2str((tot_est_cells1+tot_est_cells2)./scan_length_min),'\n'));
fprintf(strcat('\nMaximum channel peak count rate (/hr):\t\t', num2str(max(cor_pks1,cor_pks2)*60./scan_length_min),'\n'));
fprintf(strcat('Maximum channel cell detection rate (/hr):\t', num2str(max(tot_est_cells1,tot_est_cells2)*60./scan_length_min),'\n'));
fprintf(strcat('\nNumber super-peaks found:\t\t\t', num2str(num_super_peaks1+num_super_peaks2),'\n'));
fprintf(strcat('Super-peak detection rate (/min):\t\t', num2str((num_super_peaks1+num_super_peaks2)./scan_length_min),'\n'));



%% Output parameters
% Input parameters for processing
in_dat.background_sliding_window = background_sliding_window;
in_dat.smoothing_siding_window = smoothing_siding_window;
in_dat.rel_thresh = rel_thresh;
in_dat.prominence_factor = prominence_factor;
in_dat.coinc_window = coinc_window;
in_dat.est_cell_amplitude = est_cell_amplitude;
in_dat.super_peak_factor = super_peak_factor;

% Output processed data
out_dat.scan_length = time(end); % seconds
out_dat.ch1_pks = time(ch1_info.locs); % Peaks in channel 1
out_dat.ch2_pks = time(ch2_info.locs); % Peaks in channel 2
if length(out_dat.ch1_pks) >= length(out_dat.ch2_pks)
    detections = sort(out_dat.ch1_pks);
else
    detections = sort(out_dat.ch2_pks);
end
out_dat.detections = detections;


% Save to file
if save_file == 1
    pp_file = strcat(save_to, data_file_name,'_out.mat');
    save(pp_file, 'in_dat', 'out_dat');
end

end






% Sub routines (functions) used in the analysis


function [data_bs, s_dev, data_bg, ch_info]  = pre_proc_ch(ch, fs, bsw, time, rel_thresh, prominence_factor);
% This function ("pre_proc_ch") pre-process each data channel according to user-defined parameters above
% also finds peak candidates based on an adaptive threshold method
warning off
snr_win_length = 1;                                 % window length (in seconds) for sliding estimate of background standard deviation
mov_thresh_int = 60;                           % analysis interval in seconds for the adaptive threshold
sd_alpha = 0.5;                                   % smoothing factor for standard-deviation estimation

% perform basic background subtraction
data_bg = medfilt1(ch, bsw,'truncate');            % apply 'bsw' median filter,
data_bs =  ch-data_bg;                            % subtract background
% this sequence of code estimates the pre-processed signal background standard deviation (noise) over time.
% This is used to deterimine if a detected peak is above the PFA.
sig_std = movstd(data_bs,snr_win_length*fs);
std_smooth = zeros(1,length(sig_std));
std_smooth(1) =sig_std(1);
std_int = zeros(1,length(std_smooth));


% smoothing is used here because real peaks give a transient increase in standard deviation, which should not be included in estimate of noise. Included/excluded regions are stored in std_int
for indx = 2:length(std_smooth)
    if(sig_std(indx) > 1.05*(std_smooth(indx-1)))
        std_smooth(indx) = std_smooth(indx-1);
    else
        std_smooth(indx) = sd_alpha*std_smooth(indx-1)+(1-sd_alpha)*sig_std(indx);
        std_int(indx) = 100;
    end
end
std_smooth = smooth(std_smooth);
s_dev = mean(std_smooth);

% this sequence of code generates a moving detection threshold based on the
% computed standard deviation estimation and the user-defined PFA using the inverse-Q function
seg_length = mov_thresh_int*fs;         % length of segment to consider in samples
scan_length_mins = floor(length(std_smooth)/seg_length);
if(scan_length_mins >= 1)
    for k = 1:scan_length_mins
        %moving_threshold(k) = (mean(std_smooth(((k-1)*seg_length+1):k*seg_length))) *qfuncinv(pfa);
        moving_threshold(k) = (mean(std_smooth(((k-1)*seg_length+1):k*seg_length))) *rel_thresh;
    end
else
    %moving_threshold = (mean(std_smooth))*qfuncinv(pfa);
    moving_threshold = (mean(std_smooth))*rel_thresh;
end

moving_threshold = smooth(smooth(moving_threshold));
pks = [];
locs = [];
widths = [];
clear ch_info;
thresh_curve = zeros(size(data_bs));
% search for peaks using the adaptive threshold in each interval
if(scan_length_mins >= 1)
    for pk_indx = 1:scan_length_mins
        peak_th = moving_threshold(pk_indx);
        if(pk_indx ~= scan_length_mins)
            [pk, loc, w] =...
                findpeaks(data_bs((pk_indx-1)*seg_length+1:(pk_indx)*seg_length),...
                'MinPeakHeight',peak_th,'MinPeakProminence',peak_th./prominence_factor);
            thresh_curve((pk_indx-1)*seg_length+1:(pk_indx)*seg_length) = peak_th;
        else
            [pk,loc, w]=findpeaks(data_bs((pk_indx-1)*seg_length+1:end),...
                'MinPeakHeight',peak_th,'MinPeakProminence',peak_th./prominence_factor);
            thresh_curve((pk_indx-1)*seg_length+1:end) = peak_th;
        end
        %update pks and locs data
        pks = [pks; pk];
        locs = [locs;loc+(pk_indx-1)*seg_length];
        widths = [widths; w];
    end
else
    peak_th = moving_threshold;
    [pks,locs, widths]=findpeaks(data_bs(1:end),...
        'MinPeakHeight',peak_th,'MinPeakProminence',peak_th./prominence_factor);
    thresh_curve(1:end) = peak_th;
end

ch_info.pks = pks;
ch_info.locs = locs;
ch_info.widths = widths;
ch_info.thresh_curve = thresh_curve;

end

function [data] = check_noisy_pmt(data, bsw, fiber);
% this function checks the original data files, and looks for obvious variations between the two PMTs.
% note that these are 2 PMTs from the same fiber, not 2 different fibers, hence the shape of the curves should be very similar
%Since these are measuring the same thing obvious variations are possilbly a sign of a faulty PMT or system noise
data_bs(:,1) = data(:,1)-medfilt1(data(:,1), bsw, 'truncate');
data_bs(:,2) = data(:,2)-medfilt1(data(:,2), bsw, 'truncate');

% look for peaks on both PMT channels 10X the background noise, i.e. big very peaks
ch1_noise_thresh = std(data_bs(:,1))*10;
ch2_noise_thresh = std(data_bs(:,2))*10;
[pk_amps1 locs1] = findpeaks(data_bs(:,1),'MinPeakHeight', ch1_noise_thresh, 'MinPeakProminence', ch1_noise_thresh/2);
[pk_amps2 locs2] = findpeaks(data_bs(:,2),'MinPeakHeight', ch2_noise_thresh, 'MinPeakProminence', ch2_noise_thresh/2);
pks1 = length(pk_amps1);
pks2 = length(pk_amps2);

if (pks1 ~= 0 && pks2 ~= 0 &&  pks1 >10*pks2) || (pks2 == 0 && pks1 >10)
    fprintf(strcat('\tWARNING, suspected bad PMT-1 on fiber:\t', num2str(fiber)'))
    fprintf('\n\tGenerating figure of PMT data for review...\n')
    figure; subplot(2,1,1); plot(data_bs(:,1)); xlabel('Data number (samples)'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 1, Fiber', num2str(fiber)))
    subplot(2,1,2); plot(data_bs(:,2)); xlabel('Sample Number'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 2, Fiber', num2str(fiber)))
    prompt = '\n\nDo you wish to over-write PMT-1 data with PMT-2 data? Y/N [N]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
    if str == 'Y' || str =='y'
        data(:,1) = data(:,2);
    end
    
end

if (pks1 ~= 0 && pks2 ~= 0 &&  pks2 >10*pks1) || (pks1 == 0 && pks2 >10)
    fprintf(strcat('\tWARNING, suspected bad PMT-2 on fiber:\t', num2str(fiber)'))
    fprintf('\n\tGenerating figure of PMT data for review...\n')
    figure; subplot(2,1,1); plot(data_bs(:,1)); xlabel('Data number (samples)'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 1, Fiber', num2str(fiber)))
    subplot(2,1,2); plot(data_bs(:,2)); xlabel('Sample Number'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 2, Fiber', num2str(fiber)))
    prompt = '\n\nDo you wish to over-write PMT-2 data with PMT-1 data? Y/N [N]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
    if str == 'Y' || str =='y'
        data(:,2) = data(:,1);
    end
end


end

function [ch1_info_t, ch2_info_t, coinc_pk_count, multi_pk_count] = find_coinc_peaks(ch1_info, ch2_info, fs, coinc_window);
% looking for coincident peaks in both channels, which is usually
% indicative of motion artifiacts

cws = coinc_window*fs;
ch1_pks = length(ch1_info.pks);
coinc_pk_count = 0;
% make a temporary structure to hold the updated data.
ch1_info_t = ch1_info;
ch2_info_t = ch2_info;

% find any peaks that match each other between the two fibers with the
% user-specified coincidence window
for j = 1:ch1_pks
    pk_loc = ch1_info.locs(j);
    c_ind = find(ch2_info.locs>pk_loc-cws & ch2_info.locs<pk_loc+cws);
    if(~isempty(c_ind))
        % if any are found zero them out
        coinc_pk_count = coinc_pk_count +1;
        ch1_info_t.pks(j) = 0;
        ch1_info_t.locs(j) = 0;
        ch1_info_t.widths(j) = 0;
        ch2_info_t.pks(c_ind) = 0;
        ch2_info_t.locs(c_ind) = 0;
        ch2_info_t.widths(c_ind) = 0;
    end
end

% remove all the zeros from ch1_info and ch2_info
ch1_info_t.pks = ch1_info_t.pks(find(ch1_info_t.pks ~= 0));
ch1_info_t.locs = ch1_info_t.locs(find(ch1_info_t.locs ~= 0));
ch1_info_t.widths = ch1_info_t.widths(find(ch1_info_t.widths~= 0));
ch2_info_t.pks = ch2_info_t.pks(find(ch2_info_t.pks ~= 0));
ch2_info_t.locs = ch2_info_t.locs(find(ch2_info_t.locs ~= 0));
ch2_info_t.widths = ch2_info_t.widths(find(ch2_info_t.widths~= 0));
%now check channel 1 to see if there are any incidences of multiple peaks inside cws

multi_pk_count = 0;
ch1_pks = length(ch1_info_t.pks);
for j = 1:ch1_pks
    pk_loc = ch1_info_t.locs(j);
    if pk_loc ~= 0
        c_ind = find(ch1_info_t.locs>pk_loc-cws & ch1_info_t.locs<pk_loc+cws);
        % if there are any bunched up peaks, peak the largest one and zero out the others
        if(length(c_ind)>1)
            multi_pk_count = multi_pk_count +1;
            biggest_pk = find(ch1_info_t.pks(c_ind) == max(ch1_info_t.pks(c_ind)));
            biggest_pk_ind = c_ind(biggest_pk);
            if(length(biggest_pk_ind)>1); biggest_pk_ind = biggest_pk_ind(1);   end
            
            not_biggest_pk = c_ind(c_ind~=biggest_pk_ind) ;
            ch1_info_t.pks(not_biggest_pk) = 0;
            ch1_info_t.locs(not_biggest_pk) = 0;
            ch1_info_t.widths(not_biggest_pk) = 0;
            
            
        end
    end
end
% repeat for channel 2

ch2_pks = length(ch2_info_t.pks);
for j = 1:ch2_pks
    pk_loc = ch2_info_t.locs(j);
    if pk_loc ~= 0
        c_ind = find(ch2_info_t.locs>pk_loc-cws & ch2_info_t.locs<pk_loc+cws);
        % if there are any bunched up peaks, peak the largest one and zero out the others
        if(length(c_ind)>1)
            multi_pk_count = multi_pk_count +1;
            biggest_pk = find(ch2_info_t.pks(c_ind) == max(ch2_info_t.pks(c_ind)));
            biggest_pk_ind = c_ind(biggest_pk);
            
            if(length(biggest_pk_ind)>1); biggest_pk_ind = biggest_pk_ind(1);   end
            
            not_biggest_pk = c_ind(c_ind~=biggest_pk_ind) ;
            ch2_info_t.pks(not_biggest_pk) = 0;
            ch2_info_t.locs(not_biggest_pk) = 0;
            ch2_info_t.widths(not_biggest_pk) = 0;
        end
    end
end

% remove all the zeros from ch1_info and ch2_info
ch1_info_t.pks = ch1_info_t.pks(find(ch1_info_t.pks ~= 0));
ch1_info_t.locs = ch1_info_t.locs(find(ch1_info_t.locs ~= 0));
ch1_info_t.widths = ch1_info_t.widths(find(ch1_info_t.widths~= 0));
ch2_info_t.pks = ch2_info_t.pks(find(ch2_info_t.pks ~= 0));
ch2_info_t.locs = ch2_info_t.locs(find(ch2_info_t.locs ~= 0));
ch2_info_t.widths = ch2_info_t.widths(find(ch2_info_t.widths~= 0));


end

%SNR=20*log10(mean(ans.ch1_pks)/2);
