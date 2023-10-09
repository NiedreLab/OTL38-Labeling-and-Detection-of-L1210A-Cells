function [filteredPeaks] = combinePeaks(peaks1, peaks2)
%FILTERPEAKS Returns a peak struct with the peaks specified by keptPeaks
% keptPeaks is a logical array or array of indicies
filteredPeaks = peaks1;
filteredPeaks.pks = [peaks1.pks; peaks2.pks];
filteredPeaks.locs = [peaks1.locs; peaks2.locs];
filteredPeaks.widths = [peaks1.widths; peaks2.widths];
filteredPeaks.proms = [peaks1.proms; peaks2.proms];
filteredPeaks.count = length(filteredPeaks.pks);

end

