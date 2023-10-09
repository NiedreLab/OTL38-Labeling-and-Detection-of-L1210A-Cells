function [filteredPeaks] = filterPeaks(peaks, keptPeaks)
%FILTERPEAKS Returns a peak struct with the peaks specified by keptPeaks
% keptPeaks is a logical array or array of indicies
filteredPeaks = peaks;
filteredPeaks.pks = filteredPeaks.pks(keptPeaks);
filteredPeaks.locs = filteredPeaks.locs(keptPeaks);
filteredPeaks.widths = filteredPeaks.widths(keptPeaks);
filteredPeaks.proms = filteredPeaks.proms(keptPeaks);
filteredPeaks.count = length(filteredPeaks.pks);

end

