function DefectStructure = CalculateAngularMSD(DefectStructure)
% Computes the Angular Mean Squared Displacement (MSAD) using orientation
% vectors u(t) = [cos(phi), sin(phi)] for positive defects only.
%
% Inputs:
%   DefectStructure : structure containing defect trajectories with fields
%       .phi      : orientation angle of the defect over time (radians)
%       .lifetime : number of frames the defect exists
%       .charge   : defect charge
%
% Outputs:
%   DefectStructure : updated structure containing
%       .AngularMSD : angular mean squared displacement as a function of lag time
%
% Notes:
%   - The MSAD is computed using orientation vectors to avoid angular
%     discontinuities (e.g. jumps at ±π).
%   - Only positive defects are considered because negative defects do not
%     have a well-defined polarity.
%   - Lag times are computed up to half of the defect lifetime to ensure
%     sufficient statistics.

for def = 1:numel(DefectStructure)

    % ---- Process only positive defects ----
    if DefectStructure(def).charge <= 0
        DefectStructure(def).AngularMSD = [];
        continue
    end

    philist  = DefectStructure(def).phi;
    lifetime = DefectStructure(def).lifetime;

    % Maximum lag time considered
    numberOfDeltaT = floor(lifetime/2);

    % Valid indices (frames where orientation exists)
    indices = find(~isnan(philist(:)));

    if numel(indices) < 2
        DefectStructure(def).AngularMSD = [];
        continue
    end

    % Orientation angles
    phi = philist(indices);

    % Orientation unit vectors
    u = [cos(phi), sin(phi)];

    AngularMSD = zeros(numberOfDeltaT,1);

    % --- Compute MSAD ---
    for dt = 1:numberOfDeltaT

        % Orientation vectors separated by lag time dt
        u0 = u(1:end-dt,:);
        u1 = u(1+dt:end,:);

        % Vector difference
        du = u1 - u0;

        % Mean squared angular displacement
        AngularMSD(dt) = mean(sum(du.^2,2), 'omitnan');

    end

    % Store result
    DefectStructure(def).AngularMSD = AngularMSD;

end
end