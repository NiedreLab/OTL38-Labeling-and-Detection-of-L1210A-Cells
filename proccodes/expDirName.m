function [experiment] = expDirName(expDir)
%EXPDIRNAME Name of the experiment derived from it's path
%   Detailed explanation goes here
if ispc
    experiment = expDir((find(expDir == '\', 1, 'last')+1):end);
else
    experiment = expDir((find(expDir == '/', 1, 'last')+1):end);
end

end

