function [scan] = scanDirName(scanDir)
%SCANDIRNAME Name of the scan derived from it's path
%  Useful if uigetdir() is used to select the experiment. Supports Windows
%  and UNIX filepaths
if ispc
    scan = scanDir((find(scanDir == '\', 1, 'last')+1):end);
else
    scan = scanDir((find(scanDir == '/', 1, 'last')+1):end);
end

end

