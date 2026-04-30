function [ThetaValues, PhiValues, ThetaNegative, AnglesNegative] = ...
    AngleDistributionv2_withCoords(Matrix, xMat, yMat, offsetrow, offsetcol)
% Computes the angular distribution around a central point using the
% orientation field stored in Matrix and real-space coordinates.
%
% Inputs:
%   Matrix     : 2D matrix containing the orientation angle at each grid point (rad)
%   xMat       : 2D matrix containing the x-coordinate of each grid node
%   yMat       : 2D matrix containing the y-coordinate of each grid node
%   offsetrow  : row offset applied to the center position
%   offsetcol  : column offset applied to the center position
%
% Outputs:
%   ThetaValues     : angular deviation θ (deg) for all sampled points,
%                     sorted according to the corresponding geometric angle
%   PhiValues       : geometric angle φ (deg) of each sampled point relative
%                     to the center
%   ThetaNegative   : θ values corresponding only to negative φ
%   AnglesNegative  : negative φ values (deg)
%
% Notes:
%   - The center of the analysis is defined as the middle of the matrix
%     shifted by (offsetrow, offsetcol).
%   - Neighboring points are selected within annular regions of increasing
%     radius around the center using real-space distances computed from
%     xMat and yMat.
%   - The angle φ corresponds to the geometric angle between the center
%     and each neighbor.
%   - The quantity θ = acos(|sin(Matrix - φ)|) measures the angular
%     deviation between the local orientation and the radial direction.
%   - Positive and negative φ values are separated to allow asymmetric
%     angular analysis around the center.



    [nRows, nCols] = size(Matrix);

    centerRow = ceil(nRows / 2) + offsetrow;
    centerCol = ceil(nCols / 2) + offsetcol;

    centerX = xMat(centerRow, centerCol);
    centerY = yMat(centerRow, centerCol);

    ThetaValues = [];
    PhiValues = [];
    AnglesPositive = [];
    ThetaPositive = [];
    AnglesNegative = [];
    ThetaNegative = [];

    % ---- POUR QUIVER ----
    Qx = [];
    Qy = [];
    U  = [];
    V  = [];
    AngleColor = [];

    for r = 1:6
        [neighborRows, neighborCols, angles] = ...
            FindNeighborsWithAnglesCoords(centerX, centerY, r, nRows, nCols, xMat, yMat);

        theta = zeros(length(angles),1);

        for k = 1:length(angles)
            row = neighborRows(k);
            col = neighborCols(k);

            phi = angles(k);
            matrixValue = Matrix(row, col);

            theta(k) = acos(abs(sin(matrixValue - phi)));

            % ===== CHAMP DE VECTEURS =====
            Qx(end+1) = xMat(row,col);
            Qy(end+1) = yMat(row,col);

            U(end+1)  = cos(matrixValue);
            V(end+1)  = sin(matrixValue);

            AngleColor(end+1) = phi; % angle géométrique
        end

        theta  = rad2deg(theta);
        angles = rad2deg(angles);

        for i = 1:numel(angles)
            if angles(i) >= 0
                AnglesPositive(end+1) = angles(i);
                ThetaPositive(end+1)  = theta(i);
            else
                AnglesNegative(end+1) = angles(i);
                ThetaNegative(end+1)  = theta(i);
            end
        end
    end

    PhiValues   = [AnglesPositive, AnglesNegative];
    ThetaValues = [ThetaPositive, ThetaNegative];

    [PhiValues, idxSort] = sort(PhiValues);
    ThetaValues = ThetaValues(idxSort);

   
end





function [rows, cols, angles] = FindNeighborsWithAnglesCoords(centerX, centerY, radius, nRows, nCols, xMat, yMat)
    rows = [];
    cols = [];
    angles = [];

    for r = 1 : nRows
        for c = 1 : nCols
            % Real-space distance using xMat, yMat
            dx_real = xMat(r, c) - centerX;
            dy_real = yMat(r, c) - centerY;
            dist_real = sqrt(dx_real^2 + dy_real^2);

            % Keep if within annulus thickness ±1
            if dist_real < radius + 1 && dist_real > radius - 1
                rows(end+1) = r;
                cols(end+1) = c;
                theta = atan2(dy_real, dx_real);
                angles(end+1) = theta;
            end
        end
    end
end
