function pairs = pairPositiveDefects_phiproduct(defects, tol)
%PAIRPOSITIVEDEFECTS_PHIPRODUCT Pair positive defects and compute polarity projection.
%
%   pairs = pairPositiveDefects_phiproduct(defects, tol)
%
%   This function identifies all pairs of positive defects (+1/2 charge)
%   at each time point and computes their relative distance and polarity
%   projection along the connecting axis.
%
%   INPUTS
%   defects : structure array
%       Structure containing tracked defects. Each element must contain:
%           - charge     : defect topological charge
%           - positions  : (T x 2) matrix of defect coordinates [row, col]
%           - phi        : (T x 1) orientation angle of the defect (radians)
%
%   tol : double (optional)
%       Numerical tolerance used to avoid division by zero when computing
%       normalized vectors. Default value is 1e-6.
%
%   OUTPUT
%   pairs : structure array (length = number of time points)
%       For each time point t:
%           pairs(t).time  : time index
%           pairs(t).list  : structure array containing pair information
%
%       Each pair entry contains:
%           - idx1      : index of the first defect
%           - idx2      : index of the second defect
%           - distance  : distance between defects (µm)
%           - rpp       : polarity projection along the pair axis
%
%   The polarity projection rpp is defined as:
%
%       rpp = dot(r_hat , (p1 - p2)/2)
%
%   where r_hat is the normalized vector connecting the two defects,
%   and p1, p2 are the polarity vectors derived from their orientation
%   angles.

if nargin < 2
    tol = 1e-6;
end

% Determine maximum number of time points
nT = 0;
for d = 1:numel(defects)
    nT = max(nT, size(defects(d).positions,1));
end

pairs = struct('time', cell(nT,1), 'list', cell(nT,1));

for t = 1:nT

    list = struct('idx1', {}, 'idx2', {}, ...
                  'distance', {}, 'rpp', {});

    posIdxAll = find([defects.charge] > 0);

    % ---- Filter valid positive defects ----
    validPos = [];

    for k = 1:numel(posIdxAll)

        idx = posIdxAll(k);

        if size(defects(idx).positions,1) < t
            continue
        end

        pos = defects(idx).positions(t,:);
        if all(pos == 0)
            continue
        end

        phi = defects(idx).phi;
        if isempty(phi) || t > numel(phi) || isnan(phi(t))
            continue
        end

        validPos(end+1) = idx; %#ok<AGROW>
    end

    nPos = numel(validPos);

    if nPos < 2
        pairs(t).time = t;
        pairs(t).list = [];
        continue
    end

    % ---- Loop over all pairs of positive defects ----
    for i = 1:nPos
        for j = i+1:nPos

            idx1 = validPos(i);
            idx2 = validPos(j);

            P1 = defects(idx1).positions(t,:);
            P2 = defects(idx2).positions(t,:);

            r = P2 - P1;
            r_norm = norm(r);

            if r_norm < tol
                continue
            end

            r_hat = r / r_norm;

            % defect polarity vectors
            phi1 = defects(idx1).phi(t);
            phi2 = defects(idx2).phi(t);

            p1 = [cos(phi1), sin(phi1)];
            p2 = [cos(phi2), sin(phi2)];

            % polarity projection along connecting axis
            rpp = dot(r_hat, (p1 - p2)/2);

            % store pair
            list(end+1).idx1     = idx1; %#ok<AGROW>
            list(end).idx2       = idx2;
            list(end).distance   = r_norm * 30 * 0.65;
            list(end).rpp        = rpp;

        end
    end

    pairs(t).time = t;
    pairs(t).list = list;

end

end