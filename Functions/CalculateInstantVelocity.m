function DefectStructure = CalculateInstantVelocity(DefectStructure, pixelSize, nodeSpacing, dt_min)
%CALCULATEINSTANTVELOCITY Compute instantaneous and mean defect velocities.
%
%   DefectStructure = CalculateInstantVelocity(DefectStructure, pixelSize, nodeSpacing, dt_min)
%
%   This function computes the instantaneous displacement and velocity of
%   each tracked defect between consecutive frames, as well as the total
%   displacement and mean velocity over its lifetime.
%
%   INPUTS
%   DefectStructure : structure array
%       Structure containing tracked defects. Each element must include:
%           - positions : (T x 2) array of defect coordinates [row, column]
%
%   pixelSize : double
%       Size of one pixel in micrometers (µm).
%
%   nodeSpacing : double
%       Spatial scaling factor applied to convert grid spacing into pixels.
%
%   dt_min : double
%       Time interval between frames in minutes.
%
%   OUTPUT
%   DefectStructure : structure array
%       Same structure as input, with additional fields:
%           - dx_inst_um            : instantaneous x-displacement (µm)
%           - dy_inst_um            : instantaneous y-displacement (µm)
%           - v_inst_um_per_h       : instantaneous velocity magnitude (µm/h)
%           - totalDisplacement_um  : total vector displacement [dx dy] (µm)
%           - totalDistance_um      : total path length traveled (µm)
%           - meanVelocity_um_per_h : mean velocity averaged over valid steps (µm/h)
%
%   Instantaneous velocities are computed between consecutive valid frames
%   using Euclidean displacement and converted to µm/h.

    dt_hours = dt_min / 60;              % convert minutes to hours
    scale = pixelSize * nodeSpacing;     % µm per pixel

    for def = 1:numel(DefectStructure)

        positions = DefectStructure(def).positions; % [x y]
        x = positions(:,1);
        y = positions(:,2);

        N = numel(x);

        % Initialize
        dx = nan(N-1,1);
        dy = nan(N-1,1);
        validStep = false(N-1,1);

        for t = 1:(N-1)

            % Skip invalid positions
            if any(isnan([x(t), y(t), x(t+1), y(t+1)])) || ...
               (x(t)==0 && y(t)==0) || ...
               (x(t+1)==0 && y(t+1)==0)
                continue
            end

            dx(t) = x(t+1) - x(t);
            dy(t) = y(t+1) - y(t);
            validStep(t) = true;
        end

        % Convert to microns
        dx_um = dx * scale;
        dy_um = dy * scale;

        % Instantaneous displacement magnitude
        instDistance_um = sqrt(dx_um.^2 + dy_um.^2);

        % Instantaneous velocity
        v_inst_um_per_h = instDistance_um / dt_hours;

        % Total displacement vector
        total_dx = nansum(dx_um);
        total_dy = nansum(dy_um);

        % Total traveled distance
        totalDistance_um = nansum(instDistance_um);

        % Mean velocity
        meanVel = mean(v_inst_um_per_h(validStep), 'omitnan');

        % Store results
        DefectStructure(def).dx_inst_um = dx_um;
        DefectStructure(def).dy_inst_um = dy_um;
        DefectStructure(def).v_inst_um_per_h = v_inst_um_per_h;
        DefectStructure(def).totalDisplacement_um = [total_dx total_dy];
        DefectStructure(def).totalDistance_um = totalDistance_um;
        DefectStructure(def).meanVelocity_um_per_h = meanVel;

    end
end