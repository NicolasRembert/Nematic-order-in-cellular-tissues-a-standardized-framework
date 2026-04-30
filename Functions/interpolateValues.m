function Interp = interpolateValues(x, y, x_0, y_0, Values)
    % interpolateValues: Interpolates the given values using scatteredInterpolant.
    %
    % Inputs:
    % - x: Target x-coordinates for interpolation (matrix or vector)
    % - y: Target y-coordinates for interpolation (matrix or vector)
    % - x_0: Original x-coordinates of the data points
    % - y_0: Original y-coordinates of the data points
    % - Values: Values at (x_0, y_0) to be interpolated
    %
    % Output:
    % - Interp: Interpolated values at (x, y)

    % Ensure Values is a column vector
    values = double(reshape(Values, [], 1));

    % Create a scatteredInterpolant object with the 'nearest' method
    F = scatteredInterpolant(x_0, y_0, values, 'nearest');

    % Perform interpolation at the target points
    Interp = F(x, y);
end
