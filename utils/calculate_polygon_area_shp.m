function areas_m2 = calculate_polygon_area_shp(S)
%CALCULATE_POLYGON_AREA Compute polygon areas (m²) for WGS84 geographic shapefiles.
%   AREAS_M2 = CALCULATE_POLYGON_AREA(S) computes the area of each polygon
%   feature in a struct array returned by SHAPEREAD with 'UseGeoCoords', true.
%   Coordinates are assumed to be in WGS84 geographic (lat/lon in degrees).
%
%   INPUT
%       S         - Struct array from:
%                   shaperead(..., 'UseGeoCoords', true)
%                   Each element must contain either:
%                       - Lat / Lon fields (degrees), or
%                       - Y / X fields   (degrees, interpreted as lat/lon).
%
%   OUTPUT
%       areas_m2  - N×1 numeric vector of polygon areas in square meters,
%                   where N = numel(S). The i-th entry corresponds to S(i).
%
%   REQUIREMENTS
%       - Mapping Toolbox (functions: WGS84ELLIPSOID, AREAINT).
%       - Coordinates must be WGS84 geographic (latitude/longitude in degrees).
%
%   BEHAVIOUR
%       - NaN vertices (segment separators or missing data) are removed
%         before area computation.
%       - Polygons with fewer than 3 valid vertices are assigned area 0.
%
%   EXAMPLE
%       % Read polygons in geographic coordinates (lat/lon, degrees):
%       S = shaperead('admin_boundaries.shp', 'UseGeoCoords', true);
%       % Compute areas:
%       A = calculate_polygon_area(S);
%       % Convert to km²:
%       A_km2 = A / 1e6;
%
%   SEE ALSO
%       shaperead, wgs84Ellipsoid, areaint

    % ------------------------ Input checks --------------------------------
    if nargin < 1
        error('calculate_polygon_area:NotEnoughInputs', ...
            'Input struct array S is required.');
    end

    if ~isstruct(S)
        error('calculate_polygon_area:InvalidInput', ...
            'Input S must be a struct array (e.g., from SHAPEREAD).');
    end

    % Detect coordinate field names (Lat/Lon or Y/X)
    hasLatLon = isfield(S, 'Lat') && isfield(S, 'Lon');
    hasYX     = isfield(S, 'Y')   && isfield(S, 'X');

    if ~hasLatLon && ~hasYX
        error('calculate_polygon_area:MissingFields', ...
            'Input S must contain Lat/Lon or Y/X fields (in degrees).');
    end

    if hasLatLon
        latCell = {S.Lat};
        lonCell = {S.Lon};
    else
        latCell = {S.Y};
        lonCell = {S.X};
    end
    areas_m2 = calculate_polygon_area_lonslats_cell(latCell, lonCell);  
    
end
