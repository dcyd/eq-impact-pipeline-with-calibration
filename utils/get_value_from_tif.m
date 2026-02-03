function point_values = get_value_from_tif(points, value_data, R)
% GET_value_FROM_TIF Retrieve .tif values from a GeoTIFF file.
%
%   gmf_values = GET_VALUE_FROM_TIF(points, GMF_data, R)
%
%   This function assigns values to a set of points based on
%   the values defined in a GeoTIFF file representing the distribution of the value.
%
%   INPUTS:
%       points      - Nx2 array of [longitude, latitude] coordinates.
%       value_data, R - the variable after readgeoraster the path to the GeoTIFF file.
%
%   OUTPUT:
%       gmf_values  - Nx1 array of values corresponding to
%                     the input points.
%
%   EXAMPLE:
%       % Define points (longitude and latitude)
%       points = [
%           -118.25, 34.05;   % Los Angeles
%           -122.42, 37.77;   % San Francisco
%           -121.89, 37.33;   % San Jose
%           -117.16, 32.72;   % San Diego
%       ];
%
%       % Get ground motion values from the GeoTIFF file
%       [gmf_data,gmf_R] = readgeoraster('ground_motion_distribution.tif');
%       gmf_values = get_value_from_tif(points,gmf_data,gmf_R);
%
%       % Display the results
%       for i = 1:length(gmf_values)
%           fprintf('Point (%.4f, %.4f) has GMF value: %.3f\n', ...
%                   points(i,1), points(i,2), gmf_values(i));
%       end

    % Initialize the output with NaN values for unmatched points
    num_points = size(points, 1);
    point_values = NaN(num_points, 1);
    
    % Extract point coordinates
    x_points = points(:, 1);
    y_points = points(:, 2);

    % Transform points to pixel coordinates based on the referencing matrix R
    % pixel_coords = worldToPixel(R, x_points, y_points);
    [pixel_coords_x,pixel_coords_y] = geographicToIntrinsic(R,y_points,x_points);
    % Round pixel coordinates to nearest integer
    row_indices = round(pixel_coords_y);
    col_indices = round(pixel_coords_x);
    
    % Create logical indices to check if pixel coordinates are within bounds
    valid_indices = (row_indices > 0 & row_indices <= size(value_data, 1) & ...
                     col_indices > 0 & col_indices <= size(value_data, 2));
    
    % Assign GMF values only for valid indices
    point_values(valid_indices) = value_data(sub2ind(size(value_data), ...
                                                   row_indices(valid_indices), ...
                                                   col_indices(valid_indices)));
end
