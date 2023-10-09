function [updatedFwdPeaks,updatedRevPeaks, updatedFwdScores, updatedRevScores, updatedFwdSpeeds, updatedRevSpeeds, doubleMatchCount] = determineBestDirectionalMatch(fwdPeaks,revPeaks, fwdScores, revScores,fwdSpeeds, revSpeeds,probe)
% If there are peaks that are both fwd and rev matched, determineBestDirectionalMatch
% chooses best directional match based on the match scoring
%Josh Pace 20210905


%Input arguments:
%'fwdPeaks' is the forward peaks struct
%'revPeaks' is the reverse peaks struct
%'fwdScores' is the forward matched scores
%'revScores' is the reverse matched scores
%'fwdSpeeds' is the forward matched speeds
%'revSpeeds' is the reverse matched speeds
%'probe' is 1 or 2 to specify if the function looks at probe 1 or 2 
%peak locations in the 'fwdPeaks' and 'revPeaks' structs 


%Output argument
%'updatedFwdPeaks' is the forward peaks struct after checking matches
%'updatedRevPeaks' is the reverse peaks struct after checking matches
%'updatedFwdScores' is the forward matched scores after checking matches
%'updatedRevScores' is the reverse matched scores after checking matches
%'updatedFwdSpeeds' is the forward matched speeds after checking matches
%'updatedRevSpeeds' is the reverse matched speeds after checking matches

%"double matches" is when the same peaks are matched in fwd and reverse
%directions



%Getting logical arrays which will allow only lower scoring double matches
%to remain
fwdKeptPeaks = logical(ones(1,length(fwdPeaks(probe).locs)));
revKeptPeaks = logical(ones(1,length(revPeaks(probe).locs)));

%Getting array element locations of double matches by using peak
%locations
fwdDoubleMatches = ismember(fwdPeaks(probe).locs,revPeaks(probe).locs);
revDoubleMatches = ismember(revPeaks(probe).locs,fwdPeaks(probe).locs);

%For the case when there are no double matches, function will just return
%the input
if nnz(fwdDoubleMatches) == 0 || nnz(revDoubleMatches) == 0
    updatedFwdPeaks = fwdPeaks;
    updatedRevPeaks = revPeaks;
    updatedFwdScores = fwdScores;
    updatedRevScores = revScores;
    updatedFwdSpeeds = fwdSpeeds;
    updatedRevSpeeds = revSpeeds;
    doubleMatchCount = 0;
    return
end

%Pick out array indices of double matches
[fwdInds,~,~] = find(fwdDoubleMatches);
[revInds,~,~] = find(revDoubleMatches);

%Double checking that number of double matches are the same in each
%direction (should be)
if length(fwdInds) ~= length(revInds)
    disp('Index arrays not same length!')
end

%For each double match case, the fwd and rev scores are compared and the
%direction in each the score is lower is kept
%fwd and rev peak locations are double checked to verify that the correct
%postions and scores are being checked
for i = 1:length(fwdInds)
    fwdScoreValue = fwdScores(fwdInds(i));
    fwdPeakLocValue = fwdPeaks(probe).locs(fwdInds(i));
    revScoreValue = revScores(revInds(i));
    revPeakLocValue = revPeaks(probe).locs(revInds(i));
    
    if fwdPeakLocValue == revPeakLocValue
    
        if fwdScoreValue < revScoreValue
            revKeptPeaks(revInds(i)) = 0; %remove reverse match case
        else
            fwdKeptPeaks(fwdInds(i)) = 0; %remove fwd match case
        end
    else
        disp("Locations are not the same!")
    end
    
    
    
end
    %Update fwd and rev peak structs to remove 
    updatedFwdPeaks(1) = filterPeaks(fwdPeaks(1), fwdKeptPeaks);
    updatedFwdPeaks(2) = filterPeaks(fwdPeaks(2), fwdKeptPeaks);
    updatedRevPeaks(1) = filterPeaks(revPeaks(1), revKeptPeaks);
    updatedRevPeaks(2) = filterPeaks(revPeaks(2), revKeptPeaks);
    %Update fwd and rev score arrays
    updatedFwdScores = fwdScores(fwdKeptPeaks);
    updatedRevScores = revScores(revKeptPeaks);
    %Update fwd and rev speed arrays
    updatedFwdSpeeds = fwdSpeeds(fwdKeptPeaks);
    updatedRevSpeeds = revSpeeds(revKeptPeaks);
    %Getting the number of double matches that were identified
    doubleMatchCount = length(fwdInds);

end

