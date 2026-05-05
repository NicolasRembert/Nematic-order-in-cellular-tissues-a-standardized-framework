function DefectStructure = CalculateMSD(DefectStructure, pixelSize, nodeSpacing)
%CALCULATEMSD Compute mean squared displacement (MSD) of tracked defects
%
% DefectStructure = CalculateMSD(DefectStructure, pixelSize, nodeSpacing)
%
% Inputs
%   DefectStructure - structure array containing defect trajectories
%                     (must include the field .positions)
%   pixelSize       - physical size of one pixel
%   nodeSpacing     - spacing factor between nodes in the grid
%
% Outputs
%   DefectStructure - updated structure array with additional fields
%                     .MSD              : mean squared displacement curve
%                     .alpha_firstHalf  : MSD scaling exponent at short times
%                     .alpha_secondHalf : MSD scaling exponent at longer times
%
% The function computes the mean squared displacement of each defect
% trajectory and estimates the scaling exponent α from log–log fits
% of MSD vs time at early and late times.

    nIgnore = 0;

    for def = 1:numel(DefectStructure)

        positions = DefectStructure(def).positions;

        if size(positions,1) <= nIgnore + 5

            DefectStructure(def).MSD = NaN;
            DefectStructure(def).alpha_firstHalf  = NaN;
            DefectStructure(def).alpha_secondHalf = NaN;
            continue

        end

        positions = positions(nIgnore+1:end,:);

        indices = find(positions(:,1) ~= 0);
        data = positions(indices,:);

        numberOfDeltaT = floor(size(data,1)*0.5);

        if numberOfDeltaT < 5

            DefectStructure(def).MSD = NaN;
            DefectStructure(def).alpha_firstHalf  = NaN;
            DefectStructure(def).alpha_secondHalf = NaN;
            continue

        end

        MSD = nan(numberOfDeltaT,1);

        %% ================= MSD =================

        for dt = 1:numberOfDeltaT

            dx = data(1+dt:end,1) - data(1:end-dt,1);
            dy = data(1+dt:end,2) - data(1:end-dt,2);

            squaredDisplacement = ...
                sum((dx * pixelSize * nodeSpacing).^2) + ...
                sum((dy * pixelSize * nodeSpacing).^2);

            MSD(dt) = squaredDisplacement / length(dx);

        end


        %% ================= ALPHA =================

        t = (1:numel(MSD))';
        maxDt = length(MSD);

        % Adaptive fitting ranges

        earlyRange = 2:min(10,maxDt);

        if maxDt >= 20
            lateRange = 10:20;
        elseif maxDt > 10
            lateRange = 10:maxDt;
        else
            lateRange = [];
        end


        % -------- EARLY --------

        if ~isempty(earlyRange)

            x1 = log(t(earlyRange));
            y1 = log(MSD(earlyRange));

            valid1 = isfinite(x1) & isfinite(y1) & y1 > 0;

            if sum(valid1) >= 2

                p1 = polyfit(x1(valid1),y1(valid1),1);
                alpha_first = p1(1);

                if alpha_first < 0
                    alpha_first = NaN;
                end

            else
                alpha_first = NaN;
            end

        else
            alpha_first = NaN;
        end


        % -------- LATE --------

        if ~isempty(lateRange)

            x2 = log(t(lateRange));
            y2 = log(MSD(lateRange));

            valid2 = isfinite(x2) & isfinite(y2) & y2 > 0;

            if sum(valid2) >= 2

                p2 = polyfit(x2(valid2),y2(valid2),1);
                alpha_second = p2(1);

                if alpha_second < 0
                    alpha_second = NaN;
                end

            else
                alpha_second = NaN;
            end

        else
            alpha_second = NaN;
        end


        %% ================= STORE =================

        DefectStructure(def).MSD               = MSD;
        DefectStructure(def).alpha_firstHalf   = alpha_first;
        DefectStructure(def).alpha_secondHalf  = alpha_second;

    end
end