function angleBox = extractSquareBox(matrix, x, y, r)
%EXTRACTSQUAREBOX Extract a square region centered on a pixel.
%
%   angleBox = extractSquareBox(matrix, x, y, r)
%
%   This function extracts a square subregion of size (2r+1) × (2r+1)
%   centered at the coordinates (x, y) from the input matrix. If the box
%   extends outside the matrix boundaries, the missing values are filled
%   with NaN.
%
%   INPUTS
%   matrix : 2D array
%       Input matrix from which the square region is extracted.
%
%   x : integer
%       X-coordinate (column index) of the center pixel.
%
%   y : integer
%       Y-coordinate (row index) of the center pixel.
%
%   r : integer
%       Radius of the square box. The final box size is (2r+1) × (2r+1).
%
%   OUTPUT
%   angleBox : (2r+1) × (2r+1) array
%       Extracted square region centered at (x, y). Values outside the
%       original matrix boundaries are filled with NaN.

    % Get matrix size
    [x_end, y_end] = size(matrix);

    % Determine box boundaries inside the matrix
    x1 = max(1, x - r);
    x2 = min(x_end, x + r);
    y1 = max(1, y - r);
    y2 = min(y_end, y + r);

    % Initialize output box
    angleBox = NaN(2*r + 1, 2*r + 1);

    % Corresponding indices in the original matrix
    rowRange = y1:y2;
    colRange = x1:x2;

    % Indices in the output box
    rowIdxAngleBox = (y1 - y + r + 1):(y2 - y + r + 1);
    colIdxAngleBox = (x1 - x + r + 1):(x2 - x + r + 1);

    % Copy values
    angleBox(rowIdxAngleBox, colIdxAngleBox) = matrix(rowRange, colRange);

end