% Script to begin NIR DiFC data collection, saving, and plotting
%
% ------- CHANGELOG ----------------------
% V. Pera, X.Tan; June 2016

% X. Tan; Nov. 2016, add two channels

% X. Tan; May.2017, check if the file exist

% mdshumel; 22.Jan.2021, branched from red DiFC. Modified system for single
% channel detection for NIR DiFC prototype

% mdshumel; 22.Jan.2021, update organization, file naming, update comments,
% fix typos

% mdshumel; 25.Jan.2021, wrap channel name, yLabel, igain, and plot
% linespec into the channels structure

% mdshumel; 28.Jan.2021, updated to work with universal scripts

% mdshumel; 29.Jan.2021, added press enter to stop DAQ session

% mdshumel; 29.Jan.2021, forked to work with blue-green DiFC system

% ------- FILE DEPENDENCIES --------------
% DiFC_SaveData.m
% DiFC_RTMAPlot.m
% StopDiFC.m

close all; clear; clc;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be set by User: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where the saving/plotting funcs are located
daqcodespath = 'C:\Users\Mark Niedre\Northeastern University\Niedre, Mark - Niedre_Lab\DiFC Files\codes\daqcodes';

% Where data will be saved
saveDir = 'C:\Users\Mark Niedre\Northeastern University\Niedre, Mark - Niedre_Lab\Malcolm\MATLAB\codes\blue-green\test';
stem = 'bg_test_1';

% DAQ Session properties
samplingRate = 2000; % in Hz
saveInterval = 10; % how often to save data (in s)
filterTime = 0.01;   % how long in time the moving average window should be (in s)
disptime = 20; % how much time real-time plot should show (in s)

% Channels structure contains four fields:
channels.name = '';     % name of channel
channels.units = '';    % units of signal
channels.igain = [];    % scales raw data to proper unit
channels.linespec = ''; % specifies the linespec on the generated plots
channels.file = [];     % specifies which file (..._F1 - Fn_...) to save channel

channels(1).name = 'Fiber1, Collection1';
channels(1).units = 'mV';
channels(1).igain = 1000;
channels(1).linespec = '-r';
channels(1).file = 1;

channels(2).name = 'Fiber1, Collection2';
channels(2).units = 'mV';
channels(2).igain = 1000;
channels(2).linespec = '-b';
channels(2).file = 1;

channels(3).name = 'Fiber2, Collection1';
channels(3).units = 'mV';
channels(3).igain = 1000;
channels(3).linespec = '-k';
channels(3).file = 2;

channels(4).name = 'Fiber2, Collection2';
channels(4).units = 'mV';
channels(4).igain = 1000;
channels(4).linespec = '-g';
channels(4).file = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if save folder exists
if exist(saveDir, 'dir') == 0
    disp('Error: folder does not exist')
    return
end

% Check if funcpath exists
if exist(daqcodespath, 'dir') == 0
    disp('Error: could not find functions')
    return
end

% Check if the data already exists
if saveDir(end) ~= '\'; saveDir = [saveDir, '\']; end
fname = [saveDir, stem];
if exist([fname, '_F1_1.mat'],'file') ~= 0
    disp('Error: the file already exists')
    return
end

% TODO: Ensure this folder a) exists and b) has required files
% Add funcs folder to path
addpath(daqcodespath);

% Set defaults
for i = 1:length(channels)
    if isempty(channels(i).name); channels(i).name = 'Unnamed Channel'; end
    if isempty(channels(i).units); channels(i).units = 'V'; end
    if isempty(channels(i).igain); channels(i).igain = 1; end
    if isempty(channels(i).linespec); channels(i).linespec = '-b'; end
    if isempty(channels(i).file); channels(i).file = 1; end
end

disp('Acquiring data in background: Press enter to stop')

% Set up DAQ Session
s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0:length(channels)-1, 'Voltage');
s.Rate = samplingRate;
s.NotifyWhenDataAvailableExceeds = s.Rate/5; % 5 times per sec
%s.DurationInSeconds = 30; save data for 30s and stop automatically
s.IsContinuous = true;  % continue acquiring data indefinitely until stop command is issued

window = filterTime*samplingRate;

lh_s = addlistener(s,'DataAvailable', @(src, event) DiFC_SaveData(src, event, fname, saveInterval*5, window, disptime, channels));

% Start data acquisition, wait for enter, stop data acquisition
s.startBackground();
pause
StopDiFC
