function [ out_dat ] = DiFC_fastBG_matching( data_path, data_file, save_flag, save_to, plot_flag)
%DiFC_process_LLC -  Data processing routine for DiFC data, specificially
%the 'voltage output' system with LLC tumor metatasis model
% M.Niedre May 2019; Revised Oct 2019
% code written to analyze DiFC data in LLC lung metastasis model.
% the code expects two files corresponding to measurements from fiber 1 and 2, respectively (_F1, _F2).
% Revised Oct 2019: allows removal of first 10 minutes of data. 
% Allows matching of peaks between channels

%%
%Note that this code inverts the data from _F1 and _F2 files
%%   User-Defined Parameters  

plot_flag = 1;
if plot_flag == 1
    figure;
end

%  Parameters for data pre-processing
background_sliding_window = 2.5;                % length of median filter used for background subtraction in seconds
smoothing_siding_window = 0.001;                % length of median filter used for data filtering in seconds
remove_10min =  0;                                      % flag to remove first 10 minutues of data (=1) which can be noisy
limit_to_35min = 0;

%   Parameters for  peak detection routine
rel_thresh = 5;                                     % multiplicative factor for the adaptive threshold;  threshold is rel_thresh*estimated standard deviation
prominence_factor = 1;                            % minimum peak prominence (versus adjacent peaks) is the threshold/prominence_factor 
max_speed = 300;                                  % maximum speed (mm/s) permitted for spike matching routine

%% Data Pre-procesisng and Peak Candidate Detection
%   Reads in data from the two fibers, performs pre-processing, and peak  detection
fprintf('\n\n\n\n\n\n\nDiFC Data Processing Routine for LLC Mice - May 2019\n')
%loads the user-supplied raw data (data_file) file for fiber 1
dat_f   = strcat(data_path, data_file, '_F1.mat'); 
load(dat_f);
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
[data] = check_noisy_pmt(data, bsw, 1, rel_thresh);

% pre-process the sum of channel 1 and 2 for fiber 1.  Finds peak
% candidates based on a calculated adaptive threshold.
fprintf('\tPre-processing data for Fiber 1 and looking for peaks\n')
data_sum = data(:,1)+data(:,2);
if remove_10min == 1
    fprintf('\tWarning: First 10 minutes removed from analysis.  If not intended clear flag "remove_10min" \n')
    data_sum = data_sum(60*10/dt+1:end);
    time = time(1:end-60*10/dt);
end
if limit_to_35min == 1 && length(data_sum) > 35*60*fs
    fprintf('\tWarning: First at most 50 min analyzed.  If not intended clear flag "limit_to_50min" \n')
    data_sum = data_sum(1:35*60*fs);
    time = time(1:35*60*fs);
end
data_sum = smooth(data_sum, ssw);
[ch_pp(:,1), ch_sd(1), ch_bg(:,1), ch1_info]  = pre_proc_ch(data_sum, fs, bsw, time, rel_thresh, prominence_factor);
clear data;


%load and processes the raw data for fiber 2 and find peak candidates (repeat of the above steps for the second data file)
dat_f   = strcat(data_path, data_file, '_F2.mat'); 
load(dat_f);
%%%%%%%%%%%%%%
% data = -data;
%%%%%%%%%%%%%%
fprintf(strcat('\nLoading raw data file for fiber 2:\t', data_file, '\n'))
if exist('time') ~= 1 || exist('data') ~= 1
    fprintf('Error: time and data variables not found - check raw data file\n')
    return
end

%check for noisy/bad PMT measurement (occasionally happens)
fprintf('\tChecking for bad/noisy PMTs on Fiber 2\n')
[data] = check_noisy_pmt(data, bsw ,2, rel_thresh);
fprintf('\tPre-processing data for Fiber 2 and looking for peaks\n')
data_sum = data(:,1)+data(:,2);
if remove_10min == 1
    data_sum = data_sum(60*10/dt+1:end);
    time = time(1:end-60*10/dt);
end
if limit_to_35min == 1 && length(data_sum) > 35*60*fs
    data_sum = data_sum(1:35*60*fs);
    time = time(1:35*60*fs);
end
data_sum = smooth(data_sum, ssw);
[ch_pp(:,2), ch_sd(2), ch_bg(:,2), ch2_info]  = pre_proc_ch(data_sum, fs, bsw, time, rel_thresh, prominence_factor);
clear data;

fprintf('\n\n')
scan_length_min = time(end)/60;
fprintf(strcat('Scan Length in Minutes:\t\t\t\t', num2str(scan_length_min),'\n'));

%fprintf(strcat('Channel 1 Mean estimated standard deviation:\t', num2str(ch_sd(1)),'\n'));
ch1_pk_cand = length(ch1_info.locs);
fprintf(strcat('Channel 1 Peak candidates identified:\t\t', num2str(length(ch1_info.locs)),'\n'));

%fprintf(strcat('Channel 2 Mean Estimated Standard Deviation:\t', num2str(ch_sd(2)),'\n'));
ch2_pk_cand = length(ch2_info.locs);
fprintf(strcat('Channel 2 Peak candidates identified:\t\t', num2str(length(ch2_info.locs)),'\n'));

% generate a plot of the background-subtracted signal, peak candidates,
% and the calculated adaptive threshold.
    ch1_info.locs = int32(ch1_info.locs);
    ch2_info.locs = int32(ch2_info.locs);
if plot_flag == 1
    subplot(2,1,1); plot(time, ch_pp(:,1)); hold on;
    plot(time(ch1_info.locs), ch_pp(ch1_info.locs,1), 'go', 'MarkerSize', 6);
    plot(time, ch1_info.thresh_curve, 'g-'); xlabel('Time (s)'); ylabel('DiFC Signal (mV)');  title('Processed DiFC Data'); legend('DiFC Data', 'Peak Candidates', 'Dynamic Threshold')
    subplot(2,1,2); plot(time, ch_pp(:,2)); hold on;
    plot(time(ch2_info.locs), ch_pp(ch2_info.locs,2), 'go', 'MarkerSize', 6 ); plot(time, ch2_info.thresh_curve, 'g-'); xlabel('Time (s)'); ylabel('DiFC Signal (mV)'); legend('DiFC Data', 'Peak Candidates', 'Dynamic Threshold')
end

%% Part 2 - Search for forward and reversed matched peaks between the two channels ("temporal coincidence detection")

% extract some basic information about the peak candidates
ch1_pks = length(ch1_info.locs); ch2_pks = length(ch2_info.locs);  mat_length = max(ch1_pks, ch2_pks);  % how many peak candidates on each channel, and find the higher of the two
all_w(1:ch1_pks,1) = ch1_info.widths*dt/.001; all_w(1:ch2_pks,2) = ch2_info.widths*dt/.001;  % put all the peak candidate widths (in ms) into one matrix
all_pk_times(1:ch1_pks,1) = time(ch1_info.locs); all_pk_times(1:ch2_pks,2) = time(ch2_info.locs);   % put all the peak candidate times (s) into one matrix
all_pks(1:ch1_pks,1) = ch1_info.pks; all_pks(1:ch2_pks,2) = ch2_info.pks;  % put all the peak candidate amplitues into one matrix

% Call the peak-matching routine in the forward direction (channel 1-->2);
fprintf('\nMatching peaks in forward direction, Channel 1-->2\t\t\t\t\t')
[true_counts, match_data, match_ampl, match_fib_speed, match_score] = peak_match(all_w, all_pk_times, all_pks, ch_sd, mat_length, max_speed);
% gather the indices of peakes that match
mpi = find(match_data ~= 0);

% Next call peak-matching routine in the reverse direction (channel 2-->1);
fprintf('\nMatching peaks in reverse direction, Channel 2-->1\t\t\t\t\t');

% reverse the channel 1 and 2 data and recall the routine
fprintf('\n\nSearching for incidences of peak matches in both directions and if found correcting')
all_w_rev = fliplr(all_w); all_pk_times_rev = fliplr(all_pk_times); all_pks_rev = fliplr(all_pks);
[true_counts_rev, match_data_rev, match_ampl_rev, match_fib_speed_rev, match_score_rev] = peak_match(all_w_rev, all_pk_times_rev, all_pks_rev, ch_sd, mat_length, max_speed);
mpi_rev = find(match_data_rev ~= 0);
fprintf(strcat('\n\tInitial count estimate Ch 1-->2:\t', num2str(true_counts)))
fprintf(strcat('\n\tInitial count estimate Ch 2-->1:\t', num2str(true_counts_rev)))

%Look for double-matches in the forward and reverse directions; if they exist chose the one with the lower match score.
num_matches = true_counts;
for j = 1:num_matches
    cand1_t = all_pk_times(mpi(j),1);  doub_match_ind = [];
    doub_match_ind = find(match_data_rev(mpi_rev) == cand1_t);
    if ~isempty(doub_match_ind)
        score1 = match_score(mpi(j));
        score2 = match_score_rev(mpi_rev(doub_match_ind));
        if score1 < score2
            match_data_rev(mpi_rev(doub_match_ind)) = 0; mpi_rev(doub_match_ind) = 1;
            true_counts_rev = true_counts_rev - 1;
        else
            match_data(mpi(j)) = 0; mpi(j) = 1;
            true_counts = true_counts - 1;
        end
    end
end
% output corrected matched counts in forward and reverse directions
fprintf(strcat('\n\tCorrected counts Ch 1-->2:\t\t', num2str(true_counts)))
fprintf(strcat('\n\tCorrected counts Ch 2-->1:\t\t', num2str(true_counts_rev)))

% Clear zero entries from mpi and mpi_rev
mpi(find(mpi == 1)) = []; mpi_rev(find(mpi_rev == 1)) = [];
if plot_flag == 1
    if true_counts>=1
        subplot(2,1,1); hold on; plot(all_pk_times(mpi,1), all_pks(mpi,1), 'r>', 'MarkerSize', 8); legend('DiFC Data', 'Peak Candidates', 'Dynamic Threshold', 'Match 1-->2')
        subplot(2,1,2); hold on; plot(match_data(mpi), match_ampl(mpi), 'r>', 'MarkerSize', 8); legend('DiFC Data', 'Peak Candidates', 'Dynamic Threshold', 'Match 1-->2')
    end
    if true_counts_rev >=1
        subplot(2,1,1); hold on; plot(match_data_rev(mpi_rev), match_ampl_rev(mpi_rev), 'b<', 'MarkerSize', 8); legend('DiFC Data', 'Peak Candidates', 'Dynamic Threshold', 'Match 1-->2', 'Match 2-->1')
        subplot(2,1,2); hold on; plot(all_pk_times(mpi_rev,2), all_pks(mpi_rev,2), 'b<', 'MarkerSize', 8); legend('DiFC Data', 'Peak Candidates', 'Dynamic Threshold', 'Match 1-->2', 'Match 2-->1')
    end
end
%% Part 3 - Output some useful data

% Output a bunch of useful statistics
fprintf('\n\nForward and reverse matching statistics:')

if true_counts ~= 0
    fprintf(strcat('\n\tMatched peaks found Ch1-->2:\t\t\t\t\t', num2str(true_counts)))
    fprintf(strcat('\n\tCount rate per minute Ch1-->2:\t\t\t\t\t', num2str(true_counts./(scan_length_min))))
    fprintf('\n')
else
    fprintf('\n\tNo matched peaks found in forward direction')
    fprintf('\n')
    warning off
end

if true_counts_rev ~=0
    fprintf(strcat('\n\tMatched peaks found Ch2-->1:\t\t\t\t\t', num2str(true_counts_rev)))
    fprintf(strcat('\n\tCount rate per minute Ch2-->1:\t\t\t\t\t', num2str(true_counts_rev./scan_length_min)))
    fprintf('\n\n')
else
    fprintf('\n\tNo matched peaks found in reverse direction')
    fprintf('\n')
    warning off
end

% Outputs to save
out_dat.time = time;
out_dat.scan_length = scan_length_min;
out_dat.ch1_pk_cand = ch1_pk_cand;
out_dat.ch2_pk_cand = ch2_pk_cand;
out_dat.cells_arterial = all_pk_times(mpi,1);
out_dat.cells_venous = all_pk_times(mpi_rev,2);
out_dat.avg_cells_arterial = length(all_pk_times(mpi,1)) / scan_length_min;
out_dat.avg_cells_venous = length(all_pk_times(mpi_rev,2)) / scan_length_min;
out_dat.ch1_match_ampl = all_pks(mpi,1);

in_dat.rel_thresh           = rel_thresh;
in_dat.prominence_factor             = prominence_factor;
in_dat.background_sliding_window    = background_sliding_window;
in_dat.smoothing_siding_window      = smoothing_siding_window;
in_dat.max_speed =      max_speed;

if save_flag == 1
    save_file = strcat(save_to, data_file,'_out.mat');
    save(save_file, 'out_dat', 'in_dat');
end

end






% Sub routines (functions) used in the analysis


function [data_bs, s_dev, data_bg, ch_info]  = pre_proc_ch(ch, fs, bsw, time, rel_thresh, prominence_factor)
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

function [data] = check_noisy_pmt(data, bsw, fiber, rel_thresh)
% this function checks the original data files, and looks for obvious variations between the two PMTs.  
% note that these are 2 PMTs from the same fiber, not 2 different fibers, hence the shape of the curves should be very similar
%Since these are measuring the same thing obvious variations are possilbly a sign of a faulty PMT or system noise
data_bs(:,1) = data(:,1)-medfilt1(data(:,1), bsw, 'truncate');
data_bs(:,2) = data(:,2)-medfilt1(data(:,2), bsw, 'truncate');

% look for peaks on both PMT channels peak_threshX the background noise, i.e. big very peaks
% peak_thresh = 10;
peak_thresh = rel_thresh;
ch1_noise_thresh = std(data_bs(:,1))*peak_thresh;
ch2_noise_thresh = std(data_bs(:,2))*peak_thresh;
[pk_amps1 locs1] = findpeaks(data_bs(:,1),'MinPeakHeight', ch1_noise_thresh, 'MinPeakProminence', ch1_noise_thresh/2); 
[pk_amps2 locs2] = findpeaks(data_bs(:,2),'MinPeakHeight', ch2_noise_thresh, 'MinPeakProminence', ch2_noise_thresh/2); 
pks1 = length(pk_amps1);
pks2 = length(pk_amps2);

noisy_thresh = 3;
if (pks1 ~= 0 && pks2 ~= 0 &&  pks1 > noisy_thresh*pks2) || (pks2 == 0 && pks1 > noisy_thresh)
        fprintf(strcat('\tWARNING, suspected bad PMT-1 on fiber:\t', num2str(fiber)'));
        fprintf('\nCh1 peaks: %d\nCh2 peaks: %d', pks1, pks2);
        fprintf('\n\tGenerating figure of PMT data for review...\n');
        figure; subplot(2,1,1); plot(data_bs(:,1)); hold on; xl = xlim; plot(xl, ch1_noise_thresh*[1 1]); xlabel('Data number (samples)'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 1, Fiber', num2str(fiber)));
        plot(locs1, pk_amps1, 'o');
        subplot(2,1,2); plot(data_bs(:,2)); hold on; xl = xlim; plot(xl, ch2_noise_thresh*[1 1]); xlabel('Sample Number'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 2, Fiber', num2str(fiber)))
        plot(locs2, pk_amps2, 'o');
%         prompt = '\n\nDo you wish to over-write PMT-1 data with PMT-2 data? Y/N [N]: ';
%         str = input(prompt,'s');
        str = 'n';
        if isempty(str)
              str = 'N';
        end
        if str == 'Y' || str =='y'
            data(:,1) = data(:,2);
        end
        
end

if (pks1 ~= 0 && pks2 ~= 0 &&  pks2 > noisy_thresh*pks1) || (pks1 == 0 && pks2 > noisy_thresh)
        fprintf(strcat('\tWARNING, suspected bad PMT-2 on fiber:\t', num2str(fiber)'));
        fprintf('\nCh1 peaks: %d\nCh2 peaks: %d', pks1, pks2);
        fprintf('\n\tGenerating figure of PMT data for review...\n');
        figure; subplot(2,1,1); plot(data_bs(:,1)); hold on; xl = xlim; plot(xl, ch1_noise_thresh*[1 1]); xlabel('Data number (samples)'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 1, Fiber', num2str(fiber)))
        plot(locs1, pk_amps1, 'o');
        subplot(2,1,2); plot(data_bs(:,2)); hold on; xl = xlim; plot(xl, ch2_noise_thresh*[1 1]); xlabel('Sample Number'); ylabel('Amplitude (mV)'); title(strcat('Signal - PMT 2, Fiber', num2str(fiber)))
        plot(locs2, pk_amps2, 'o');
%         prompt = '\n\nDo you wish to over-write PMT-2 data with PMT-1 data? Y/N [N]: ';
%         str = input(prompt,'s');
        str = 'n';
        if isempty(str)
              str = 'N';
        end
        if str == 'Y' || str =='y'
            data(:,2) = data(:,1);
        end
end


end


function [true_counts, match_data, match_ampl, match_fib_speed, match_score] = peak_match(all_w, all_pk_times, all_pks, ch_sd, mat_length, max_speed)

%normalize the peak heights according to the relative channel brightness
all_pks_sc(:,2) = all_pks(:,2)*mean(all_pks(:,1))/mean(all_pks(:,2));
mst = 3/max_speed;


% look for coincidental peak in second channel (indicative of a motion
% and/or breating artifact) and zero out peak candidates
for j = 1:mat_length
     pk_cand_time = all_pk_times(j,1);
     coinc_ind = find(all_pk_times(:,2) >= pk_cand_time-mst & all_pk_times(:,2) <= pk_cand_time+mst);
     if ~isempty(coinc_ind)
         all_pks(j,1) = 0; all_w(j,1) = 0; all_pk_times(j,1) = 0;
         all_pks(coinc_ind,2) = 0; all_w(coinc_ind,2) = 0; all_pk_times(coinc_ind,2) = 0;
         coinc_ind = [];
     end 
end
    
%initialize variables (fill matrices with zeros)
match_data = zeros(mat_length,1); match_width = zeros(mat_length,1); match_height = zeros(mat_length,1); match_speed = zeros(mat_length,1); match_snr = zeros(mat_length,1); match_score = zeros(mat_length, 1); match_fib_speed = zeros(mat_length,1); match_ampl = zeros(mat_length,1);
for j = 1:mat_length
    % based on the width of the spikes in channel 1, estimate the time delay to arrive at channel 2, 3mm distance separation.
    est_dt_pk2 = 3*all_w(j,1)/1000;
    % look for corresponding peak(s) in second channel from 0.1 s to up to 3X the estimated time delay
    time_range_low = mst + all_pk_times(j,1);     time_range_high = 10*est_dt_pk2 + all_pk_times(j,1);
    %all_pk_times(j,1)
    cand_pk_inds = find(all_pk_times(:,2) > time_range_low & all_pk_times(:,2) < time_range_high);
    % if there are no matches, true_counts does not change
    if isempty(cand_pk_inds)
        match_data(j) = 0;
    % if there is one or more matches, consider the similarity of the peaks with respect to width and height
    % choose the closer match inside the range of 0.5-2
    else
        % note: abs(log-base-2) < 1 means it's within a factor of 2 in either direction
        width_ratios = log2(all_w(cand_pk_inds,2)./all_w(j,1));
        height_ratios = log2(all_pks_sc(cand_pk_inds,2)./all_pks(j,1));
        speeds_d = 3./(all_pk_times(cand_pk_inds,2) - all_pk_times(j,1));
        speeds_w = 1000./((all_w(cand_pk_inds,2)+ all_w(j,1))/2);
        speed_ratios = log2(speeds_d./speeds_w);
        score = abs(width_ratios) + abs(height_ratios) + abs(speed_ratios);
        closest_match_ind = find(score == min(score));
        closest_match = cand_pk_inds(closest_match_ind);
        avg_snr = (all_pks(j,1)/ch_sd(1)+all_pks(closest_match,2)/ch_sd(2))/2;
        %if abs(width_ratios(closest_match_ind)) <= 1 && abs(height_ratios(closest_match_ind)) <=1
        
        if score(closest_match_ind) < 4 && avg_snr >= 3 
            %all_w(closest_match,2) = NaN; 
            match_data(j) = all_pk_times(closest_match,2);
            match_width(j) = width_ratios(closest_match_ind);
            match_ampl(j) = all_pks(closest_match, 2);
            match_height(j) = height_ratios(closest_match_ind);
            match_score(j) = score(closest_match_ind);
            match_speed(j) = speed_ratios(closest_match_ind);
            match_fib_speed(j) = speeds_d(closest_match_ind);
            match_snr(j) = (all_pks(j,1)/ch_sd(1)+all_pks(closest_match,2)/ch_sd(2))/2;
        else
            match_data(j) = 0;
        end
    end
end
% check for any duplicate entries in the matched peaks
u1 = unique(match_data); u2 = histc( match_data, u1); 
duplicates = u1(find(u2>1));
%if there are duplicates, find the closer match and remove the others
if length(duplicates) >0
    % go through each duplicate
    for j = 1:length(duplicates)
        if duplicates(j) >0
         dup_ind = find(match_data == duplicates(j));                       %find the indices of the duplicate
         best_dup_ind = find(match_score == min(match_score(dup_ind)));     %find the indice of the duplicate that has the lower score (better match)
         bdi = dup_ind(find(dup_ind ~= best_dup_ind));                      %find the indeces of the duplicates the were NOT the better match
         match_data(bdi) = 0; match_width(bdi) = 0; match_height(bdi) = 0;  %clear the data 
         match_score(bdi) = 0; match_speed(bdi) = 0; match_snr(bdi) = 0;
        end
    end
end
%fprintf('\n F1 Pks\tF2 Pks\t F1-F2 Match\t Width score\t Height score\t Speed Score\t Total score\t Avg SNR \n')
%[all_pk_times(:,1) match_data match_width match_height match_speed match_score match_snr]
 

true_counts = sum(match_data>0);
end


