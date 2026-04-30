function DefectStructure = ComputeAvgSameChargeDistance(DefectStructure, PositionMatrices)
%COMPUTEAVGSAMECHARGEDISTANCE Compute distance to nearest same-charge defect.
%
%   DefectStructure = ComputeAvgSameChargeDistance(DefectStructure, PositionMatrices)
%
%   This function computes, for each defect and each time point of its
%   trajectory, the minimum distance to the closest defect with the same
%   topological charge. The computed distances are stored in the field
%   'AvgSameChargeDist' of the DefectStructure.
%
%   INPUTS
%   DefectStructure : structure array
%       Structure containing tracked defects. Each element must contain:
%           - positions : (T x 2) array of defect positions [row, column]
%           - charge    : topological charge of the defect (±0.5)
%           - lifetime  : number of time points where the defect exists
%
%   PositionMatrices : cell array
%       Cell array containing matrices of defect charges at each time
%       point. Each matrix stores the spatial distribution of defects,
%       where entries correspond to defect charges (e.g. ±0.5).
%
%   OUTPUT
%   DefectStructure : structure array
%       Same structure as input, with an additional field:
%           - AvgSameChargeDist : (T x 1) vector containing the distance to
%             the nearest defect with the same charge at each time point.
%
%   The distance is computed in pixel units using the Euclidean metric.

    % Compute the minimum distance to the closest defect of the same charge
    % for each defect and each time point
    
    numTimetot = length(PositionMatrices);
    PositionMatrices = PositionMatrices';

    for def = 1:numel(DefectStructure)
        positions = DefectStructure(def).positions;
        numTimePoints = DefectStructure(def).lifetime;

        % Valid indices (where defect exists)
        validIdx = find(~isnan(positions(:, 1)) & positions(:,1) ~= 0);

        % Initialize output for this defect
        sameChargeDistances = NaN(numTimetot, 1);

        % Loop over time points
        for t = 1:numTimePoints
            x = positions(validIdx(t), 1);   % row
            y = positions(validIdx(t), 2);   % column
            currentCharge = DefectStructure(def).charge;

            % Find same-charge defects
            indices = find(PositionMatrices{validIdx(t)} == currentCharge);

            if isempty(indices)
                sameChargeDistances(validIdx(t)) = NaN;
                continue
            end

            % Convert indices to coordinates
            [all_x, all_y] = ind2sub(size(PositionMatrices{validIdx(t)}), indices);

            % Compute distances
            dx = all_x - x;
            dy = all_y - y;
            distances = sqrt(dx.^2 + dy.^2);

            % Remove self-distance
            distances(dx==0 & dy==0) = NaN;

            % Store minimum distance
            sameChargeDistances(validIdx(t)) = min(distances(:), [], 'omitnan');
        end

        % Save in structure
        DefectStructure(def).AvgSameChargeDist = sameChargeDistances;
    end
end