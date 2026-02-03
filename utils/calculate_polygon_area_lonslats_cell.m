function areas_m2 = calculate_polygon_area_lonslats_cell(latCell, lonCell)
%CALCULATE_POLYGON_AREA_LONSLATS_CELL Polygon areas (m²) from lat/lon cell arrays (WGS84).
%   AREAS_M2 = CALCULATE_POLYGON_AREA_LONSLATS_CELL(latCell, lonCell)
%   computes the area of each polygon given its latitude/longitude
%   coordinates stored in cell arrays. Coordinates are assumed to be in
%   WGS84 geographic (degrees).
%
%   INPUTS
%       latCell   - 1×N or N×1 cell array. latCell{i} is a numeric vector
%                   of latitudes (degrees) for the i-th polygon. May
%                   contain NaNs as segment separators.
%
%       lonCell   - 1×N or N×1 cell array. lonCell{i} is a numeric vector
%                   of longitudes (degrees) for the i-th polygon. Must have
%                   the same shape as latCell and matching lengths in each
%                   cell. May contain NaNs as segment separators.
%
%   OUTPUT
%       areas_m2  - N×1 numeric vector of polygon areas in square meters.
%                   The i-th entry corresponds to latCell{i}/lonCell{i}.
%
%   REQUIREMENTS
%       - Mapping Toolbox (functions: WGS84ELLIPSOID, AREAINT).
%       - Coordinates must be WGS84 geographic (latitude/longitude, degrees).
%
%   BEHAVIOUR
%       - NaN vertices (segment separators or missing data) are removed
%         before area computation.
%       - Polygons with fewer than 3 valid vertices are assigned area 0.
%
%   EXAMPLE
%       % Suppose latCell/lonCell come from {S.Lat} / {S.Lon}:
%       S = shaperead('admin_boundaries.shp', 'UseGeoCoords', true);
%       latCell = {S.Lat};
%       lonCell = {S.Lon};
%       A = calculate_polygon_area_lonslats_cell(latCell, lonCell);
%       A_km2 = A / 1e6;
%
%   SEE ALSO
%       shaperead, wgs84Ellipsoid, areaint

    % ------------------------ Input checks --------------------------------
    if nargin < 2
        error('calculate_polygon_area_lonslats_cell:NotEnoughInputs', ...
            'Both latCell and lonCell are required.');
    end

    if ~iscell(latCell) || ~iscell(lonCell)
        error('calculate_polygon_area_lonslats_cell:InvalidInputType', ...
            'latCell and lonCell must be cell arrays.');
    end

    if numel(latCell) ~= numel(lonCell)
        error('calculate_polygon_area_lonslats_cell:SizeMismatch', ...
            'latCell and lonCell must have the same number of elements.');
    end

    n = numel(latCell);
    areas_m2 = zeros(n, 1);

    % WGS84 reference ellipsoid, units: meters
    ell = wgs84Ellipsoid('meters');

    % ------------------------ Loop over polygons --------------------------
    for i = 1:n
        lat = latCell{i}(:);
        lon = lonCell{i}(:);

        if ~isnumeric(lat) || ~isnumeric(lon)
            error('calculate_polygon_area_lonslats_cell:CellContentType', ...
                'latCell{%d} and lonCell{%d} must be numeric vectors.', i, i);
        end

        % Remove NaNs (segment separators / invalid points)
        idx = ~isnan(lat) & ~isnan(lon);
        lat = lat(idx);
        lon = lon(idx);

        if numel(lat) < 3
            areas_m2(i) = 0;
            continue;
        end

        % AREAINT: lat, lon in degrees; ellipsoid in meters.
        % Returns area in units of ell^2 → m².
        A = areaint(lat, lon, ell);  % or areaint(lat, lon, ell, 'degrees');
        areas_m2(i) = A;
    end
end
