function [distMatPos, magMatPos, distMatNeg, magMatNeg] = compute_defect_radial_magnitude( ...
    MagnitudeData, posPOSx_new, posPOSy_new, posNEGx_new, posNEGy_new, SizeImagePixels, pixelSize)
% Computes radial magnitude distributions for each defect.
% Each defect becomes a row in the output matrices.
%
% Inputs:
%   MagnitudeData : matrix containing the magnitude field (e.g. velocity magnitude)
%   posPOSx_new, posPOSy_new : positions of positive defects
%   posNEGx_new, posNEGy_new : positions of negative defects
%   SizeImagePixels : size of the original image in pixels
%   pixelSize : size of one pixel in microns (µm)
%
% Note:
%   The "_new" suffix means that defect positions were interpolated
%   onto the velocity vector grid. These interpolated positions must
%   be used to obtain correct distances relative to MagnitudeData.
%
% Outputs:
%   distMatPos : radial distances for positive defects
%   magMatPos  : magnitude values corresponding to those distances
%   distMatNeg : radial distances for negative defects
%   magMatNeg  : magnitude values corresponding to those distances

[nRows, nCols] = size(MagnitudeData);

% Conversion factor between vector grid coordinates and image pixels
step = SizeImagePixels / nRows; % distance in image pixels per grid step

% --- Preallocate lists of rows ---
distMatPos = [];
magMatPos  = [];
distMatNeg = [];
magMatNeg  = [];

% --- Positive defects ---
for d = 1:numel(posPOSx_new)

    % Position of the positive defect (interpolated on velocity grid)
    x0 = posPOSx_new(d);
    y0 = posPOSy_new(d);

    dist_row = zeros(1, nRows*nCols);
    mag_row  = zeros(1, nRows*nCols);
    idx = 1;

    for i = 1:nRows
        for j = 1:nCols

            % Radial distance from grid point to defect (µm)
            D = sqrt(((j - x0) * step)^2 + ((i - y0) * step)^2) * pixelSize;

            % Magnitude value at this grid point
            M = MagnitudeData(i,j);

            dist_row(idx) = D;
            mag_row(idx)  = M;
            idx = idx + 1;

        end
    end

    % Add this defect as a new row
    distMatPos = [distMatPos; dist_row];
    magMatPos  = [magMatPos;  mag_row];

end

% --- Negative defects ---
for d = 1:numel(posNEGx_new)

    % Position of the negative defect (interpolated on velocity grid)
    x0 = posNEGx_new(d);
    y0 = posNEGy_new(d);

    dist_row = zeros(1, nRows*nCols);
    mag_row  = zeros(1, nRows*nCols);
    idx = 1;

    for i = 1:nRows
        for j = 1:nCols

            % Radial distance from grid point to defect (µm)
            D = sqrt(((j - x0) * step)^2 + ((i - y0) * step)^2) * pixelSize;

            % Magnitude value at this grid point
            M = MagnitudeData(i,j);

            dist_row(idx) = D;
            mag_row(idx)  = M;
            idx = idx + 1;

        end
    end

    % Add this defect as a new row
    distMatNeg = [distMatNeg; dist_row];
    magMatNeg  = [magMatNeg;  mag_row];

end

end