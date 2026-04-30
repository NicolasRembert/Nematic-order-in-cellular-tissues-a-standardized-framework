function [U_rot, V_rot, X_rot, Y_rot] = rotateMatrix(U, V, angle_rad)
%ROTATE_MATRIX Rotate a vector field and its coordinate grid around the center.
%
%   [U_rot, V_rot, X_rot, Y_rot] = rotate_matrix(U, V, angle_rad)
%
%   INPUTS:
%     U, V      : vector field components (same size)
%     angle_rad : rotation angle in radians (counterclockwise)
%
%   OUTPUTS:
%     U_rot, V_rot : rotated vector components
%     X_rot, Y_rot : rotated coordinate grids (for plotting)
%
%   The rotation is done around the matrix center, preserving the
%   coordinate structure (no image interpolation).

    if nargin < 3
        error('Usage: rotate_matrix(U, V, angle_rad)');
    end

    [rows, cols] = size(U);

    % --- Create coordinate grid ---
    [X, Y] = meshgrid(1:cols, 1:rows);
    center = [(cols+1)/2, (rows+1)/2];

    % --- Rotation matrix ---
    R = [cos(angle_rad), -sin(angle_rad);
         sin(angle_rad),  cos(angle_rad)];

    % --- Rotate coordinates ---
    pts = [X(:), Y(:)];
    pts_centered = pts - center;
    pts_rot = (R * pts_centered')' + center;
    X_rot = reshape(pts_rot(:,1), size(X));
    Y_rot = reshape(pts_rot(:,2), size(Y));

    % --- Rotate vectors ---
    uv = [U(:), V(:)];
    uv_rot = (R * uv')';
    U_rot = reshape(uv_rot(:,1), size(U));
    V_rot = reshape(uv_rot(:,2), size(V));

end


