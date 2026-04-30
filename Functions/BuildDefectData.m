function DefectData = BuildDefectData(Qmat_list, Phi_list, maxRadius, maxGap)
%BUILDDEFECTDATA Track ±1/2 topological defects across time
%
% DefectData = BuildDefectData(Qmat_list, Phi_list, maxRadius, maxGap)
%
% Inputs
%   Qmat_list   - cell array of charge matrices (values typically -0.5, 0, 0.5)
%   Phi_list    - cell array of orientation matrices corresponding to each frame
%   maxRadius   - maximum distance allowed for defect matching between frames
%   maxGap      - maximum number of consecutive missing frames before a defect track is terminated
%
% Outputs
%   DefectData  - structure containing tracked defects and associated data
%                 .defects       : structure array with defect trajectories
%                 .Qpos          : cell array of +1/2 defect masks
%                 .Qneg          : cell array of -1/2 defect masks
%                 .numTimePoints : number of frames
%                 .maxRadius     : matching distance used for tracking
%                 .maxGap        : allowed gap length in tracking
%
% The function separates positive and negative defects and performs
% robust tracking across frames using a distance-based assignment.

numTimePoints = length(Qmat_list);

%% 1) Strict separation ±0.5
Qpos = cell(numTimePoints,1);
Qneg = cell(numTimePoints,1);

for t = 1:numTimePoints
    mat = Qmat_list{t};
    Qpos{t} = (mat == 0.5);
    Qneg{t} = (mat == -0.5);
end

%% 2) Robust tracking (charges separated)
defects = struct('name',{},'charge',{},'positions',{},...
                 'phi',{},'missCount',{},'alive',{});

defects = TrackCharge_WithMemory(defects, Qpos, 0.5, maxRadius, Phi_list, maxGap);
defects = TrackCharge_WithMemory(defects, Qneg, -0.5, maxRadius, Phi_list, maxGap);

%% 3) Final structure
DefectData = struct();
DefectData.defects = defects;
DefectData.Qpos = Qpos;
DefectData.Qneg = Qneg;
DefectData.numTimePoints = numTimePoints;
DefectData.maxRadius = maxRadius;
DefectData.maxGap = maxGap;

end


function defects = TrackCharge_WithMemory(defects, Qcells, charge, maxRadius, Phi_Matrices, maxGap)
%TRACKCHARGE_WITHMEMORY Track defects of a given charge across frames
%
% defects = TrackCharge_WithMemory(defects, Qcells, charge, maxRadius, ...
%                                  Phi_Matrices, maxGap)
%
% Inputs
%   defects       - structure array containing existing defect tracks
%   Qcells        - cell array of logical matrices indicating defect positions
%   charge        - defect charge to track (+0.5 or -0.5)
%   maxRadius     - maximum allowed displacement between frames
%   Phi_Matrices  - cell array of orientation matrices
%   maxGap        - maximum number of missing frames allowed
%
% Outputs
%   defects       - updated structure array with tracked defects
%
% The function performs frame-to-frame tracking using a distance-based
% cost matrix and Hungarian assignment.

numTimePoints = length(Qcells);
bigCost = 1e6;

for t = 1:numTimePoints
    
    % 1) Retrieve active defects of the same charge only
    prevPos = [];
    activeIdx = [];
    
    for d = 1:length(defects)
        if defects(d).alive && defects(d).charge == charge
            
            if t == 1
                lastPos = defects(d).positions(1,:);
            else
                lastPos = defects(d).positions(t-1,:);
            end
            
            if any(lastPos ~= 0)
                prevPos = [prevPos; lastPos];
                activeIdx = [activeIdx; d];
            end
        end
    end
    
    % 2) Candidates at frame t
    [rows, cols] = find(Qcells{t});
    candidates = [rows, cols];
    
    % Case 1: no active tracks -> create all
    if isempty(activeIdx)
        for i = 1:size(candidates,1)
            newID = length(defects) + 1;
            
            defects(newID).name = sprintf('defect %d', newID);
            defects(newID).charge = charge;
            defects(newID).positions = zeros(numTimePoints,2);
            defects(newID).phi = NaN(numTimePoints,1);
            defects(newID).missCount = 0;
            defects(newID).alive = true;
            
            defects(newID).positions(t,:) = candidates(i,:);
            defects(newID).phi(t) = ...
                Phi_Matrices{t}(candidates(i,1), candidates(i,2));
        end
        continue
    end
    
    % Case 2: no candidates -> increase miss count
    if isempty(candidates)
        for k = 1:length(activeIdx)
            defectID = activeIdx(k);
            defects(defectID).missCount = defects(defectID).missCount + 1;
            
            if defects(defectID).missCount > maxGap
                defects(defectID).alive = false;
            end
        end
        continue
    end
    
    % 3) Cost matrix
    costMatrix = pdist2(prevPos, candidates);
    costMatrix(costMatrix > maxRadius) = bigCost;
    
    % 4) Global assignment (Hungarian)
    [assignments, ~] = matchpairs(costMatrix, bigCost);
    
    assignedTracks = false(length(activeIdx),1);
    assignedCandidates = false(size(candidates,1),1);
    
    % Update assigned tracks
    for k = 1:size(assignments,1)
        
        r = assignments(k,1);
        c = assignments(k,2);
        
        if costMatrix(r,c) >= bigCost
            continue
        end
        
        defectID = activeIdx(r);
        
        defects(defectID).positions(t,:) = candidates(c,:);
        defects(defectID).phi(t) = ...
            Phi_Matrices{t}(candidates(c,1), candidates(c,2));
        
        defects(defectID).missCount = 0;
        
        assignedTracks(r) = true;
        assignedCandidates(c) = true;
    end
    
    % Unassigned tracks -> miss
    for k = 1:length(activeIdx)
        if ~assignedTracks(k)
            defectID = activeIdx(k);
            defects(defectID).missCount = defects(defectID).missCount + 1;
            
            if defects(defectID).missCount > maxGap
                defects(defectID).alive = false;
            end
        end
    end
    
    % Unassigned candidates -> new defects
    for c = 1:size(candidates,1)
        if ~assignedCandidates(c)
            
            newID = length(defects) + 1;
            
            defects(newID).name = sprintf('defect %d', newID);
            defects(newID).charge = charge;
            defects(newID).positions = zeros(numTimePoints,2);
            defects(newID).phi = NaN(numTimePoints,1);
            defects(newID).missCount = 0;
            defects(newID).alive = true;
            
            defects(newID).positions(t,:) = candidates(c,:);
            defects(newID).phi(t) = ...
                Phi_Matrices{t}(candidates(c,1), candidates(c,2));
        end
    end
    
end

end