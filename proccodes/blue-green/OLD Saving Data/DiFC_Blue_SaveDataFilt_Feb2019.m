
%
% Script to save data in chunks in background (4 channels) 
% and plot moving-average filtered data in real time.
%
% V. Pera,X.Tan; June 15, 2016
% X. Tan; Nov. 2016,  add two channels; May 5,2017, check if the file exist
% X. Tan; Feb. 2019, edit for use on voltage PMTs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be set by User: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dir = 'C:\Users\Mark Niedre\Documents\MATLAB\Data_for_2DiFC\Lin_lab\190830_prx1Mice_Day4\Mouse1_blood\'; % directory where data files will be saved
% dir = 'C:\Users\Mark Niedre\Documents\MATLAB\Data _for_2DiFC\Lin_Lab\190516_Phantom\';
% dir = 'C:\Users\Mark Niedre\Documents\MATLAB\Data_for_2DiFC\Niedre\CartTest_2019_10_31\';
 dir = 'C:\Users\Mark Niedre\Northeastern University\Niedre, Mark - Niedre_Lab\Fernando\CURRENT EXPERIMENTS\DG5\0.3mm SD-with_FC_validation\0.75mm deep\trial3\';
%dir = 'C:\Users\Mark Niedre\Documents\MATLAB\AmberTest\';
%stem = 'test2';

% stem = 'BT_Mouse2Blood_1mL_p55V_N3'; % file name
% stem = 'Ph_p75mm_p5V_DG4_N2'; 
 stem = 'trial3';
samplingRate = 2000; % in s
saveInterval = 10; % how often to save data (in s)
params.igain1 = 1000; % Pre-amp gain settings (Gain times mV)
params.igain2 = 1000; % Gain times mV

params.igain3 = 1000;
params.igain4 = 1000;

%params = []; % saves & plots data in Volts, without applying any gains
filterTime = 0.1;   % how long in time the moving average window should be (in s)
disptime = 20; % how much time real-time plot should show (in s)
channels(1).name = 'Fiber1, Collection1'; % struct of data names for display on plots
channels(2).name = 'Fiber1, Collection2';

channels(3).name = 'Fiber2, Collection1';
channels(4).name = 'Fiber2, Collection2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first check if the data alreasy exist 
if exist([dir,stem,'_F1_1.mat'],'file') ~= 0 
     fprintf('Error: the file is already exist\n')
     return
 end

disp('Acquiring data in background; type stopDiFC_Blue_Feb2019 to stop.')

s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0:3, 'Voltage');

s.Rate = samplingRate;

s.NotifyWhenDataAvailableExceeds = s.Rate/5; % 5 times per sec

%s.DurationInSeconds = 30; save data for 30s and stop automatically
s.IsContinuous = true;  % continue acquiring data indefinitely until stop command is issued

fname = [dir,stem];

window = filterTime*samplingRate;

lh_s = addlistener(s,'DataAvailable', @(src, event) DiFC_Blue_ABsaveDataFilt_Feb2019(src, event, fname, saveInterval*5, params, window, disptime, channels));

s.startBackground();
