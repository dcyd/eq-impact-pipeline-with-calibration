function BuildingInfo_quick = generate_building_quick_sim_data(BuildingInfo, SeismicScenario, types, isGpu)
%GENERATE_BUILDING_QUICK_SIM_DATA Precompute compact inputs for fast building simulations
%   BuildingInfo_quick = GENERATE_BUILDING_QUICK_SIM_DATA(BuildingInfo, SeismicScenario)
%   reorganizes ground‐shaking and ground‐failure fields so that buildings
%   sharing identical intensity and fragility parameters point to a common
%   simulation record, reducing memory and computation cost.
%
%   BuildingInfo_quick = GENERATE_BUILDING_QUICK_SIM_DATA(..., types)
%   restricts which groups of precomputed values are created. types is a
%   cell array of strings chosen from:
%       'GSSimValues'   - Ground shaking IM + fragility combinations
%       'GFLatSimValues' - Lateral spreading PGD + liquefaction probability
%       'GFSetSimValues' - Ground settlement PGD + liquefaction probability
%       'GFLandSimValues' - Landslide PGD + landslide probability
%   By default, all four types are generated.
%
%   BuildingInfo_quick = GENERATE_BUILDING_QUICK_SIM_DATA(..., types, isGpu)
%   sets isGpu = true to treat SeismicScenario and fragility parameters as
%   gpuArray objects (when available). If false or omitted, CPU arrays are
%   used.
%
%   INPUTS
%       BuildingInfo     - Structure with building_data and fra_params
%                          (fra_inten_code, fra_media, fra_disper,
%                          rate_cd2collap, rate_collap2fatali)
%       SeismicScenario  - Structure with ground‐shaking and ground‐failure
%                          fields (PGA, SA0p3, SA0p6, SA1p0, latPGD,
%                          setPGD, landPGD, liquProb, landProb)
%       types            - (optional) cell array of output blocks to build
%       isGpu            - (optional) logical flag to enable GPU workflow
%
%   OUTPUTS
%       BuildingInfo_quick - Structure containing one or more of:
%                            .GSSimValues
%                            .GFLatSimValues
%                            .GFSetSimValues
%                            .GFLandSimValues
%                            Each sub-structure stores de-duplicated
%                            simulation values and a bid_sid index that
%                            maps each building to a row in the unique
%                            value table.
%
%   EXAMPLE
%       % Prepare all quick-simulation blocks on CPU:
%       BIq = generate_building_quick_sim_data(BuildingInfo, SeismicScenario);
%
%       % Prepare only shaking and liquefaction-related blocks on GPU:
%       t = {'GSSimValues','GFLatSimValues'};
%       BIq = generate_building_quick_sim_data(BuildingInfo, SeismicScenario, t, true);
%
%   SEE ALSO
%       generate_road_quick_sim_data
%
%   Author: Chongyang Du
%   Email:  dcy0910@gmail.com
%   Date:   20xx-xx-xx

    % ---- Input parsing ----
    if nargin <= 2 || isempty(types)
        types = {'GSSimValues', 'GFLatSimValues', 'GFSetSimValues', 'GFLandSimValues'};
    end
    if nargin <= 3
        isGpu = false;
    end

    % Validate requested types
    valid_types = {'GSSimValues', 'GFLatSimValues', 'GFSetSimValues', 'GFLandSimValues'};
    if ~all(ismember(types, valid_types))
        error('Invalid types specified.');
    end

    % ---- Optional GPU preparation ----
    if isGpu
        fields = fieldnames(SeismicScenario);
        for i = 1:numel(fields)
            if ~isa(SeismicScenario.(fields{i}), 'gpuArray')
                SeismicScenario.(fields{i}) = single(gpuArray(SeismicScenario.(fields{i})));
            end
        end

        if ~isa(BuildingInfo.fra_params.fra_media, 'gpuArray')
            BuildingInfo.fra_params.fra_inten_code = single(gpuArray(BuildingInfo.fra_params.fra_inten_code));
            BuildingInfo.fra_params.fra_media      = single(gpuArray(BuildingInfo.fra_params.fra_media));
            BuildingInfo.fra_params.fra_disper     = single(gpuArray(BuildingInfo.fra_params.fra_disper));
        end
    end

    % ---- Unpack seismic scenario fields ----
    PGA      = SeismicScenario.PGA;
    SA0p3    = SeismicScenario.SA0p3;
    SA0p6    = SeismicScenario.SA0p6;
    SA1p0    = SeismicScenario.SA1p0;
    landPGD  = SeismicScenario.landPGD;
    latPGD   = SeismicScenario.latPGD;
    setPGD   = SeismicScenario.setPGD;
    liquProb = SeismicScenario.liquProb;
    landProb = SeismicScenario.landProb;
    clear SeismicScenario

    % ---- Basic sizes and indices ----
    N_B  = size(BuildingInfo.building_data, 1);
    N_sce = size(landPGD, 2) - 3;
    N_sta = size(BuildingInfo.fra_params.fra_media, 2); % Number of damage states

    [~, gmf_zids] = ismember(BuildingInfo.building_zone(:, 2), PGA(:, 1));
    idx = 1:N_sce;

    [zone_ids, ~, zid2build] = unique(gmf_zids);
    % zone_ids: M×1, unique zone indices
    % zid2build: N_B×1, maps each building to 1..M

    % =====================================================================
    % Ground-failure: lateral spreading
    % =====================================================================
    if ismember('GFLatSimValues', types)
        try
            [uni_values, ~, zone_sids] = unique( ...
                single([latPGD(zone_ids, idx + 3), liquProb(zone_ids, idx + 3)]), ...
                'rows');
            clear latPGD
            GFLatSimValues.latPGD   = uni_values(:, 1:N_sce);
            GFLatSimValues.liquProb = uni_values(:, N_sce+1:end);

            % For each building: get its zone, then the index in the unique table
            GFLatSimValues.bid_sid = [BuildingInfo.building_data, ...
                                      gmf_zids, zone_sids(zid2build)];
            clear uni_values
        catch ME %#ok<NASGU>
            if isa(latPGD, 'gpuArray')
                latPGD   = gather(latPGD);
                liquProb = gather(liquProb);
            end
            [uni_values, ~, zone_sids] = unique( ...
                single([latPGD(zone_ids, idx + 3), liquProb(zone_ids, idx + 3)]), ...
                'rows');
            clear latPGD
            GFLatSimValues.latPGD   = uni_values(:, 1:N_sce);
            GFLatSimValues.liquProb = uni_values(:, N_sce+1:end);
            GFLatSimValues.bid_sid  = [BuildingInfo.building_data, ...
                                       gmf_zids, zone_sids(zid2build)];
            clear uni_values
        end
        BuildingInfo_quick.GFLatSimValues = GFLatSimValues;
        clear GFLatSimValues
    end

    % =====================================================================
    % Ground-failure: settlement
    % =====================================================================
    if ismember('GFSetSimValues', types)
        try
            [uni_values, ~, zone_sids] = unique( ...
                single([setPGD(zone_ids, idx + 3), liquProb(zone_ids, idx + 3)]), ...
                'rows');
            clear setPGD liquProb
            GFSetSimValues.setPGD   = uni_values(:, 1:N_sce);
            GFSetSimValues.liquProb = uni_values(:, N_sce+1:end);
            GFSetSimValues.bid_sid  = [BuildingInfo.building_data, ...
                                       gmf_zids, zone_sids(zid2build)];
            clear uni_values
        catch ME %#ok<NASGU>
            if isa(setPGD, 'gpuArray')
                setPGD   = gather(setPGD);
                liquProb = gather(liquProb);
            end
            [uni_values, ~, zone_sids] = unique( ...
                single([setPGD(zone_ids, idx + 3), liquProb(zone_ids, idx + 3)]), ...
                'rows');
            clear setPGD liquProb
            GFSetSimValues.setPGD   = uni_values(:, 1:N_sce);
            GFSetSimValues.liquProb = uni_values(:, N_sce+1:end);
            GFSetSimValues.bid_sid  = [BuildingInfo.building_data, ...
                                       gmf_zids, zone_sids(zid2build)];
            clear uni_values
        end
        BuildingInfo_quick.GFSetSimValues = GFSetSimValues;
        clear GFSetSimValues
    end

    % =====================================================================
    % Ground-failure: landslides
    % =====================================================================
    if ismember('GFLandSimValues', types)
        try
            [uni_values, ~, zone_sids] = unique( ...
                single([landPGD(zone_ids, idx + 3), landProb(zone_ids, idx + 3)]), ...
                'rows');
            clear landPGD landProb
            GFLandSimValues.landPGD  = uni_values(:, 1:N_sce);
            GFLandSimValues.landProb = uni_values(:, N_sce+1:end);
            GFLandSimValues.bid_sid  = [BuildingInfo.building_data, ...
                                        gmf_zids, zone_sids(zid2build)];
            clear uni_values
        catch ME %#ok<NASGU>
            if isa(landPGD, 'gpuArray')
                landPGD  = gather(landPGD);
                landProb = gather(landProb);
            end
            [uni_values, ~, zone_sids] = unique( ...
                single([landPGD(zone_ids, idx + 3), landProb(zone_ids, idx + 3)]), ...
                'rows');
            clear landPGD landProb
            GFLandSimValues.landPGD  = uni_values(:, 1:N_sce);
            GFLandSimValues.landProb = uni_values(:, N_sce+1:end);
            GFLandSimValues.bid_sid  = [BuildingInfo.building_data, ...
                                        gmf_zids, zone_sids(zid2build)];
            clear uni_values
        end
        BuildingInfo_quick.GFLandSimValues = GFLandSimValues;
        clear GFLandSimValues
    end

    % =====================================================================
    % Ground shaking + fragility combinations
    % =====================================================================
    if ismember('GSSimValues', types)

        % Fragility parameters
        fra_inten_code      = BuildingInfo.fra_params.fra_inten_code;  % IM code
        fra_media           = BuildingInfo.fra_params.fra_media;       % Median
        fra_disper          = BuildingInfo.fra_params.fra_disper;      % Dispersion
        rate_cd2collap      = BuildingInfo.fra_params.rate_cd2collap;
        rate_collap2fatali  = BuildingInfo.fra_params.rate_collap2fatali;

        % Unique combinations of zone + fragility parameters
        [comb_ids, ~, bld2comb] = unique( ...
            [gmf_zids, fra_inten_code, fra_media, fra_disper, ...
             rate_cd2collap, rate_collap2fatali], ...
            'rows', 'stable');
        zone_ids  = comb_ids(:, 1); % M2×1
        int_codes = comb_ids(:, 2); % M2×1
        M2        = size(comb_ids, 1);

        try
            if isa(PGA, 'gpuArray')
                gmf_block = gpuArray.zeros(M2, N_sce, 'single');
            else
                gmf_block = zeros(M2, N_sce, 'single');
            end

            % Fill ground-motion field according to intensity measure code
            sel = (int_codes == 1);
            gmf_block(sel, :) = PGA(zone_ids(sel), idx + 3);
            clear PGA
            sel = (int_codes == 2);
            gmf_block(sel, :) = SA0p3(zone_ids(sel), idx + 3);
            clear SA0p3
            sel = (int_codes == 3);
            gmf_block(sel, :) = SA0p6(zone_ids(sel), idx + 3);
            clear SA0p6
            sel = (int_codes == 4);
            gmf_block(sel, :) = SA1p0(zone_ids(sel), idx + 3);
            clear SA1p0

            % Second row-wise unique on ground-motion + fragility combination
            [uni_vals, ~, comb2uni] = unique( ...
                single([gmf_block, comb_ids(:, 3:end)]), ...
                'rows');
            clear gmf_block

            % Assemble final output
            GSSimValues.gmf_value          = uni_vals(:, 1:N_sce);
            GSSimValues.fra_media          = uni_vals(:, N_sce + (1:N_sta));
            GSSimValues.fra_disper         = uni_vals(:, N_sce + N_sta + 1);
            GSSimValues.rate_cd2collap     = uni_vals(:, N_sce + N_sta + 2);
            GSSimValues.rate_collap2fatali = uni_vals(:, N_sce + N_sta + 3);

            % build_sids: N_B×1, maps each building to final unique row
            build_sids = comb2uni(bld2comb);
            GSSimValues.bid_sid = [BuildingInfo.building_data, build_sids];
            clear gmf_block uni_vals

        catch ME %#ok<NASGU>
            % Fallback to CPU arrays
            PGA       = gather(PGA);
            SA0p3     = gather(SA0p3);
            SA0p6     = gather(SA0p6);
            SA1p0     = gather(SA1p0);
            fra_media = gather(fra_media);
            fra_disper = gather(fra_disper);

            gmf_block = zeros(M2, N_sce, 'single');

            % Fill ground-motion field according to intensity measure code
            sel = (int_codes == 1);
            gmf_block(sel, :) = PGA(zone_ids(sel), idx + 3);
            clear PGA
            sel = (int_codes == 2);
            gmf_block(sel, :) = SA0p3(zone_ids(sel), idx + 3);
            clear SA0p3
            sel = (int_codes == 3);
            gmf_block(sel, :) = SA0p6(zone_ids(sel), idx + 3);
            clear SA0p6
            sel = (int_codes == 4);
            gmf_block(sel, :) = SA1p0(zone_ids(sel), idx + 3);
            clear SA1p0

            % Second row-wise unique on ground-motion + fragility combination
            [uni_vals, ~, comb2uni] = unique( ...
                [gmf_block, comb_ids(:, 3:end)], ...
                'rows');
            clear gmf_block

            % Assemble final output
            GSSimValues.gmf_value          = uni_vals(:, 1:N_sce);
            GSSimValues.fra_media          = uni_vals(:, N_sce + (1:N_sta));
            GSSimValues.fra_disper         = uni_vals(:, N_sce + N_sta + 1);
            GSSimValues.rate_cd2collap     = uni_vals(:, N_sce + N_sta + 2);
            GSSimValues.rate_collap2fatali = uni_vals(:, N_sce + N_sta + 3);

            % build_sids: N_B×1, maps each building to final unique row
            build_sids = comb2uni(bld2comb);
            GSSimValues.bid_sid = [BuildingInfo.building_data, build_sids];
            clear gmf_block uni_vals
        end

        BuildingInfo_quick.GSSimValues = GSSimValues;
        clear GSSimValues
    end
end
