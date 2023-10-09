function [data, time, params] = formatData(scanDir)
%FORMATDATA Takes in a scan folder and outputs data that can be immediately pre-processed
% If there are multiple PMTs on the same fiber, the data for both PMTs is
% summed

% different Oses have different slashes for file paths
if ispc
    slash = '\';
else
    slash = '/';
end

% Find all stitched files in the scan folder
files = dir([scanDir  slash '*.mat']);

pat = regexpPattern('.*(F\d\.mat)');
good_files = false(length(files), 1);
for ii = 1:length(files)
    good_files(ii) = contains(files(ii).name, pat);
end
files = files(good_files);

% If no files are found, return
if isempty(files)
    disp('Could not find stitched files')
    return
end

% Load the first file to initialize return values
load([files(1).folder slash files(1).name], 'data', 'time', 'params')

% Initialize the data/ params we intend to return
data_sum = zeros(length(time), length(files));
new_params = repmat(params(1), length(files),1);

% While we have the first file loaded, set the first column of data and
% update the first params stuct
data_sum(:,1) = sum(data,2);
new_params(1).name = sprintf('%s Fiber %g', scanDirName(scanDir), params(1).fiber);

% Update every other column of data / params struct
for ii = 2:length(files)
    load([files(ii).folder slash files(ii).name], 'data', 'time', 'params')
    data_sum(:,ii) = sum(data,2);
    new_params(ii) = params(1);
    new_params(ii).name = sprintf('%s Fiber %g', scanDirName(scanDir), params(1).fiber);
end

data = data_sum;
params = new_params;

