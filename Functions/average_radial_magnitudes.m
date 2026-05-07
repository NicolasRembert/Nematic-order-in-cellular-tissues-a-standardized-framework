function [avgDistPos, avgMagPos, stdMagPos, avgDistNeg, avgMagNeg, stdMagNeg] = average_radial_magnitudes(distMatPos, magMatPos, distMatNeg, magMatNeg)
% Averages and computes standard deviations of magnitudes for points at the same distance using nanmean and nanstd.
%
% Outputs:
%   avgDistPos - unique positive distances
%   avgMagPos  - mean magnitude at each positive distance
%   stdMagPos  - standard deviation of magnitudes at each positive distance
%   avgDistNeg - unique negative distances
%   avgMagNeg  - mean magnitude at each negative distance
%   stdMagNeg  - standard deviation of magnitudes at each negative distance

% Flatten matrices into vectors
distPosVec = distMatPos(:);
magPosVec  = magMatPos(:);
distNegVec = distMatNeg(:);
magNegVec  = magMatNeg(:);

% --- Positive defects ---
[uniqueDistPos, ~, idxPos] = unique(distPosVec);
avgMagPos = nan(size(uniqueDistPos));
stdMagPos = nan(size(uniqueDistPos));
for k = 1:numel(uniqueDistPos)
    vals = magPosVec(idxPos == k);
    avgMagPos(k) = nanmean(vals);
    stdMagPos(k) = nanstd(vals);
end
avgDistPos = uniqueDistPos;

% --- Negative defects ---
[uniqueDistNeg, ~, idxNeg] = unique(distNegVec);
avgMagNeg = nan(size(uniqueDistNeg));
stdMagNeg = nan(size(uniqueDistNeg));
for k = 1:numel(uniqueDistNeg)
    vals = magNegVec(idxNeg == k);
    avgMagNeg(k) = nanmean(vals);
    stdMagNeg(k) = nanstd(vals);
end
avgDistNeg = uniqueDistNeg;

end

