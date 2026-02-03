function [gid,gidDesc] = map_building_gid(lonLat, gidShpFile)
%MAP_BUILDING_GID Map buildings to admin units (GID_1, GID_2, GID_3).
%
%   [gid,gidDesc] = MAP_BUILDING_GID(lonLat, gidShpFile)
%
%   INPUTS
%       lonLat       - [N×2] array of [lon_deg, lat_deg]
%       gidShpFile    - path to boundary shapefile that contains
%                       fields 'GID_1', 'GID_2', 'GID_3'.
%
%   OUTPUT
%       gid     - [N×3 single], numeric IDs for each level.
%                              Column 1: GID_1 numeric code
%                              Column 2: GID_2 numeric code
%                              Column 3: GID_3 numeric code
%                              Unmapped / missing -> -1
%       gidDesc - {'gid1','gid2','gid3'}

    % if ~isfield(buildingCore, 'lonLat')
    %     error('map_building_gid:MissingLonLat', ...
    %         'buildingCore.lonLat is required.');
    % end

    coords = lonLat;          % [N×2], [lon, lat]
    N = size(coords, 1);

    % Initialize all GID levels to -1 (unmapped)
    gid_num = -ones(N, 3, 'single');       % [N×3], columns for GID_1/2/3

    % Read admin boundary shapefile
    boundaryShp = shaperead(gidShpFile, 'UseGeoCoords', true);

    % Field names for the three levels
    gidFields = {'GID_1','GID_2','GID_3'};

    % Map each GID level
    for k = 1:3
        fieldName = gidFields{k};

        % If this field does not exist in the shapefile, keep column as -1
        if ~isfield(boundaryShp, fieldName)
            continue;
        end

        % Raw values at each building coordinate (N×1 cell)
        gid_value = get_value_from_shp(coords, boundaryShp, fieldName);

        % --- Pre-classify entries: empty / numeric / string-like ---------
        isEmpty   = cellfun(@isempty, gid_value);
        isNumeric = cellfun(@isnumeric, gid_value);
        isStrLike = ~isEmpty & ~isNumeric;

        % Pre-extract numeric values (first element if vector)
        numericVals = nan(N, 1);
        if any(isNumeric&~isEmpty)
            numericVals(isNumeric&~isEmpty) = cellfun(@(v) double(v(1)), gid_value(isNumeric&~isEmpty));
        end

        % Precompute tokens for string-like entries (only once per entry)
        tokensCell = cell(N, 1);
        idxStr = find(isStrLike);
        for j = 1:numel(idxStr)
            s = char(gid_value{idxStr(j)});
            tokensCell{idxStr(j)} = regexp(s, '\d+', 'match');  % all integer substrings
        end

        % Initialize this column with -1
        col = -ones(N, 1, 'single');

        % 1) Fill numeric entries
        validNumMask = isNumeric & ~isnan(numericVals);
        if any(validNumMask)
            col(validNumMask) = single(numericVals(validNumMask));
        end

        % 2) Fill string-like entries using precomputed tokens
        for j = 1:numel(idxStr)
            idx = idxStr(j);
            tokens = tokensCell{idx};
            if isempty(tokens)
                continue;
            end

            % For GID_1, GID_2, GID_3, we typically want the 1st, 2nd, 3rd integer.
            if numel(tokens) >= k
                v = str2double(tokens{k});
            else
                % Fallback: use the last available integer
                v = str2double(tokens{end});
            end

            if ~isnan(v)
                col(idx) = single(v);
            end
        end

        gid_num(:, k) = col;
    end

    gid     = gid_num;                 % [N×3]
    gidDesc = gidFields;

end