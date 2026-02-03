function point_value = get_value_from_shp(points, polygons, fieldName)
%GET_VALUE_FROM_SHP Assign polygon attribute values to points.
%   point_value = GET_VALUE_FROM_SHP(points, polygons, fieldName)
%   returns, for each input point, the value of a specified attribute
%   from the polygon that contains that point (if any).
%
%   INPUTS
%       points    - N×2 numeric matrix of point coordinates:
%                   [longitude, latitude] (degrees) in the same coordinate
%                   system as the polygon shapefile.
%
%       polygons  - Struct array from SHAPEREAD containing polygon data.
%                   Each element is expected to have fields:
%                       .X, .Y     (polygon vertices; may include NaNs)
%                     and the attribute field specified by fieldName.
%
%       fieldName - Char or string scalar specifying the polygon attribute
%                   field to extract, e.g. 'GID_1'.
%
%   OUTPUT
%       point_value - N×1 cell array. point_value{i} is:
%                       - the attribute value of the polygon that contains
%                         points(i,:), or
%                       - [] if the point does not fall inside any polygon
%                         (or no candidate polygon passes the coarse bbox
%                         filter).
%
%   BEHAVIOUR
%       - A coarse bounding-box intersection is used first to discard
%         polygons that are far away from all points.
%       - For remaining polygons, INPOLYGON is applied to all points.
%       - If multiple polygons overlap and contain the same point, the
%         attribute of the last matching polygon in the filtered list
%         overwrites previous ones.
%
%   EXAMPLE
%       S   = shaperead('gid1_boundaries.shp', 'UseGeoCoords', true);
%       xy  = [building_data(:,3), building_data(:,4)];  % [lon,lat]
%       gid = get_value_from_shp(xy, S, 'GID_1');
%
%   SEE ALSO
%       shaperead, inpolygon

    % Extract the field values from the polygons
    fieldValues = {polygons.(fieldName)};


    % Calculate the bounding box of all points
    minLon = min(points(:,1));
    maxLon = max(points(:,1));
    minLat = min(points(:,2));
    maxLat = max(points(:,2));

    hasLatLon = isfield(polygons, 'Lat') && isfield(polygons, 'Lon');
    hasYX     = isfield(polygons, 'Y')   && isfield(polygons, 'X');

    if ~hasLatLon && ~hasYX
        error('add_building_attributes:MissingFields', ...
            'Input shapefile must contain Lat/Lon or Y/X fields (in degrees).');
    end

    if hasLatLon
        latCell = {polygons.Lat};
        lonCell = {polygons.Lon};
    else
        latCell = {polygons.Y};
        lonCell = {polygons.X};
    end

    % Extract bounding boxes for all polygons
    polyMinLon = cellfun(@min, lonCell);
    polyMaxLon = cellfun(@max, lonCell);
    polyMinLat = cellfun(@min, latCell);
    polyMaxLat = cellfun(@max, latCell);

    % Determine which polygons intersect with the points' bounding box
    bboxIntersect = (polyMaxLon >= minLon) & (polyMinLon <= maxLon) & ...
                    (polyMaxLat >= minLat) & (polyMinLat <= maxLat);

    % Filter polygons based on bounding box intersection
    polygons = polygons(bboxIntersect);
    fieldValues = fieldValues(bboxIntersect);
    lonCell = lonCell(bboxIntersect);
    latCell = latCell(bboxIntersect);

    % Initialize the output values array
    point_value = cell(size(points, 1), 1);


    % Use inpolygon for the selected polygons
    for i = 1:length(polygons)
        inPoly = inpolygon(points(:,1), points(:,2), lonCell{i}, latCell{i});
        point_value(inPoly) = fieldValues(i);
    end

end