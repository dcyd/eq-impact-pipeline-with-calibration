function [buildingGSfp, buildingGFfp] = calculate_building_damage_prob(BuildingInfo_quick, FreeParams)
%CALCULATE_BUILDING_DAMAGE_PROB Building damage exceedance probabilities.
%   [buildingGSfp, buildingGFfp] = CALCULATE_BUILDING_DAMAGE_PROB( ...
%       BuildingInfo_quick, FreeParams) computes the probabilities of
%   buildings exceeding damage states due to ground shaking and ground
%   failure under multiple seismic scenarios.
%
%   INPUTS
%       BuildingInfo_quick - Structure containing pre-processed simulation
%           inputs for each building:
%             .GSSimValues    - Ground-shaking IMs and fragility parameters
%                               (fields: gmfValue, fraMedia, fraDisper,
%                               bid_sid).
%             .GFLatSimValues - Lateral spreading PGD and liquefaction
%                               probability (fields: latPGD, liquProb,
%                               bid_sid).
%             .GFSetSimValues - Settlement PGD and liquefaction probability
%                               (fields: setPGD, liquProb, bid_sid).
%             .GFLandSimValues- Landslide PGD and landslide probability
%                               (fields: landPGD, landProb, bid_sid).
%
%       FreeParams         - (optional) structure or numeric vector of
%           correction factors applied to fragility medians:
%               If numeric (length 4):
%                   FreeParams(1) -> cfBuildingLatMedia
%                   FreeParams(2) -> cfBuildingSetMedia
%                   FreeParams(3) -> cfBuildingLandMedia
%                   FreeParams(4) -> cfBuildGSMedia
%               If structure, the following fields are recognised
%               (defaults in parentheses):
%                   .cfBuildingLatMedia  (1)
%                   .cfBuildingSetMedia  (1)
%                   .cfBuildingLandMedia (1)
%                   .cfBuildGSMedia      (1)
%                   .state               - optional string selector:
%                                          "comp": use complete damage
%                                          "exte": use extensive damage
%
%   OUTPUTS
%       buildingGSfp - Structure with ground-shaking exceedance
%           probabilities aggregated at building level:
%               .unique_GSfp - unique rows of N_B×N_sce×Ns probabilities
%               .gs_bid      - N_B×1 index mapping each building to a row
%                              in unique_GSfp.
%
%       buildingGFfp - Structure with ground-failure exceedance
%           probabilities:
%               If FreeParams.state is missing or "comp":
%                   buildingGFfp.zone_GFfp_cp - zone-level complete-damage
%                                              probabilities
%                   buildingGFfp.bid_zid      - mapping from building to
%                                              zone index
%               If FreeParams.state is "exte":
%                   buildingGFfp.zone_GFfp_ex - zone-level extensive-damage
%                                              probabilities
%                   buildingGFfp.bid_zid      - mapping from building to
%                                              zone index
%
%   EXAMPLE
%       % Default correction factors, return complete-damage GF probabilities
%       [GS, GF] = calculate_building_damage_prob(BuildingInfo_quick);
%
%       % Use custom correction factors and request extensive-damage GF
%       % probabilities
%       fp = struct('cfBuildGSMedia',1.1, ...
%                   'cfBuildingLatMedia',0.9, ...
%                   'cfBuildingSetMedia',1.0, ...
%                   'cfBuildingLandMedia',1.2, ...
%                   'state',"exte");
%       [GS, GF] = calculate_building_damage_prob(BuildingInfo_quick, fp);
%
%   SEE ALSO
%       generate_building_quick_sim_data

    % -------------------- Free parameter handling ------------------------
    if nargin < 2
        FreeParams = struct;
    end

    % Allow passing correction factors as a numeric vector
    if isnumeric(FreeParams)
        buff = FreeParams;
        FreeParams = struct;
        FreeParams.cfBuildingLatMedia  = buff(1);
        FreeParams.cfBuildingSetMedia  = buff(2);
        FreeParams.cfBuildingLandMedia = buff(3);
        FreeParams.cfBuildGSMedia      = buff(4);
    end

    % Fill in defaults if fields are missing
    if ~isfield(FreeParams, 'cfBuildingLatMedia')
        FreeParams.cfBuildingLatMedia = 1;
    end
    if ~isfield(FreeParams, 'cfBuildingSetMedia')
        FreeParams.cfBuildingSetMedia = 1;
    end
    if ~isfield(FreeParams, 'cfBuildingLandMedia')
        FreeParams.cfBuildingLandMedia = 1;
    end
    if ~isfield(FreeParams, 'cfBuildGSMedia')
        FreeParams.cfBuildGSMedia = 1;
    end

    % -------------------- Load pre-processed inputs ----------------------
    GSSimValues    = BuildingInfo_quick.GSSimValues;    % Ground shaking
    GFLatSimValues = BuildingInfo_quick.GFLatSimValues; % Lateral spreading
    GFSetSimValues = BuildingInfo_quick.GFSetSimValues; % Settlement
    GFLandSimValues= BuildingInfo_quick.GFLandSimValues;% Landslide
    clear BuildingInfo_quick

    % -------------------- Ground shaking probabilities -------------------
    buildingGSfp = calculate_building_GSfp(GSSimValues, FreeParams);

    % -------------------- Ground failure probabilities -------------------
    if ~isfield(FreeParams, 'state')
        % Default: return complete-damage probabilities
        [~, buildingGFfp] = calculate_building_GFfp( ...
            GFLatSimValues, GFSetSimValues, GFLandSimValues, FreeParams);
    else
        if strcmp(FreeParams.state, "comp")
            [~, buildingGFfp] = calculate_building_GFfp( ...
                GFLatSimValues, GFSetSimValues, GFLandSimValues, FreeParams);
        elseif strcmp(FreeParams.state, "exte")
            [buildingGFfp, ~] = calculate_building_GFfp( ...
                GFLatSimValues, GFSetSimValues, GFLandSimValues, FreeParams);
        else
            buildingGFfp = zeros(0, 0);
        end
    end
end


% ========================================================================
% Local helper: lognormal exceedance probabilities
% ========================================================================
function prob = calculateFP(gmfValue, fraMedia, fraDisper)
%CALCULATEFP Lognormal CDF of intensity vs fragility parameters.
%   prob = CALCULATEFP(gmfValue, fraMedia, fraDisper) returns the
%   probability of exceeding damage states for a given ground-motion field
%   and lognormal fragility parameters.

    is_batch   = 0;
    batch_size = 800; % Adjust based on memory limits

    % Size bookkeeping
    [N, N_sce] = size(gmfValue);
    Ns = size(fraMedia, 2);

    % Expand parameters to full size if needed
    if size(fraMedia, 1) == 1
        fraMedia = repmat(fraMedia, N, 1);
    end
    if isscalar(fraDisper)
        fraDisper = repmat(fraDisper, N, Ns);
    elseif size(fraDisper, 1) == 1
        fraDisper = repmat(fraDisper, N, 1);
    elseif size(fraDisper, 2) == 1
        fraDisper = repmat(fraDisper, 1, Ns);
    end

    isgpu = isa(gmfValue, 'gpuArray');
    if isgpu
        if ~isa(fraMedia, 'gpuArray')
            fraMedia  = single(gpuArray(fraMedia));
            fraDisper = single(gpuArray(fraDisper));
        end
    end

    fraMediaExp  = reshape(fraMedia,  N, 1, Ns); % N×1×Ns
    fraDisperExp = reshape(fraDisper, N, 1, Ns); % N×1×Ns
    clear fraMedia fraDisper

    logFraMediaExp = log(fraMediaExp);
    clear fraMediaExp

    threshold = logFraMediaExp - 3 .* fraDisperExp;
    mask      = sum(log(gmfValue) > threshold, 2) > 0;
    clear threshold

    if isgpu
        prob = gpuArray.zeros(N, N_sce, Ns, 'single');
    else
        prob = zeros(N, N_sce, Ns);
    end

    if is_batch
        num_batches = ceil(N_sce / batch_size);
        for b = 1:num_batches
            idx_start = (b - 1) * batch_size + 1;
            idx_end   = min(b * batch_size, N_sce);
            idx       = idx_start:idx_end;

            prob(mask, idx, :) = logncdf( ...
                gmfValue(mask, idx), ...
                logFraMediaExp(mask, :), ...
                fraDisperExp(mask, :));
        end
    else
        prob(mask, :, :) = logncdf( ...
            gmfValue(mask, :), ...
            logFraMediaExp(mask, :), ...
            fraDisperExp(mask, :));
    end

    clear fraDisperExp mask gmfValue logFraMediaExp
end


% ========================================================================
% Local helper: building-level ground-shaking probabilities
% ========================================================================
function buildingGSfp = calculate_building_GSfp(GSSimValues, freeparams)
%CALCULATE_BUILDING_GSFP Building-level exceedance probabilities (shaking).
%   buildingGSfp = CALCULATE_BUILDING_GSFP(GSSimValues, freeparams) returns
%   a structure with unique probability rows and a mapping from buildings
%   to these rows.

    gmfValue  = GSSimValues.gmf_value;
    fraMedia  = GSSimValues.fra_media;   % Median
    fraDisper = GSSimValues.fra_disper;  % Dispersion

    if isfield(freeparams, 'cfBuildGSMedia')
        k_bf = freeparams.cfBuildGSMedia;
    else
        k_bf = 1;
    end

    % Adjust median values
    fraMedia = fraMedia * k_bf;

    % Compute exceedance probabilities at unique IM–fragility combinations
    buildingGSfpSel = calculateFP(gmfValue, fraMedia, fraDisper);
    buildingGSfpSel = gather(buildingGSfpSel);
    GSSimValues.bid_sid = gather(GSSimValues.bid_sid);

    % Map back to building level
    buildingGSfpFull = buildingGSfpSel(GSSimValues.bid_sid(:, end), :, :);
    buildingGSfpFull(buildingGSfpFull < 1e-6) = single(0);

    [uni_building_GSfp, ~, gs_bid] = unique(buildingGSfpFull, "rows");
    clear buildingGSfpFull

    buildingGSfp.unique_GSfp = uni_building_GSfp;
    buildingGSfp.gs_bid      = gs_bid;
end


% ========================================================================
% Local helper: building-level ground-failure probabilities
% ========================================================================
function [buildingGFfpEx, buildingGFfpCp] = calculate_building_GFfp( ...
    GFLatSimValues, GFSetSimValues, GFLandSimValues, params)
%CALCULATE_BUILDING_GFFP Building-level exceedance probabilities (failure).
%   [buildingGFfpEx, buildingGFfpCp] = CALCULATE_BUILDING_GFFP(...) returns
%   structures with zone-level extensive and complete damage probabilities
%   and mappings from buildings to zones.

    if nargin < 4
        params = struct;
    end

    if isfield(params, 'cfBuildingLatMedia')
        cfBuildingLatMedia = params.cfBuildingLatMedia;
    else
        cfBuildingLatMedia = 1;
    end
    if isfield(params, 'cfBuildingSetMedia')
        cfBuildingSetMedia = params.cfBuildingSetMedia;
    else
        cfBuildingSetMedia = 1;
    end
    if isfield(params, 'cfBuildingLandMedia')
        cfBuildingLandMedia = params.cfBuildingLandMedia;
    else
        cfBuildingLandMedia = 1;
    end

    % Lateral spreading parameters
    fra_media_lat  = cfBuildingLatMedia  * 60;
    fra_disper_lat = 1.2;

    % Settlement parameters
    fra_media_set  = cfBuildingSetMedia  * 10;
    fra_disper_set = 1.2;

    % Landslide parameters
    fra_media_land  = cfBuildingLandMedia * 10;
    fra_disper_land = 0.5;

    % -------------------- Lateral failure -------------------------------
    [N, M] = size(GFLatSimValues.latPGD);
    if isa(GFLatSimValues.latPGD, 'gpuArray')
        prob_lat_fp = gpuArray.zeros(N, M, 'single');
    else
        prob_lat_fp = zeros(N, M);
    end
    mask_liq = sum(GFLatSimValues.liquProb > 0.01, 2) > 0;
    prob_lat_fp(mask_liq, :) = calculateFP( ...
        GFLatSimValues.latPGD(mask_liq, :), ...
        fra_media_lat, fra_disper_lat);
    prob_lat_fp = prob_lat_fp .* GFLatSimValues.liquProb;

    % -------------------- Settlement failure ----------------------------
    [N, M] = size(GFSetSimValues.setPGD);
    if isa(GFSetSimValues.setPGD, 'gpuArray')
        prob_set_fp = gpuArray.zeros(N, M, 'single');
    else
        prob_set_fp = zeros(N, M);
    end
    mask_liq = sum(GFSetSimValues.liquProb > 0.01, 2) > 0;
    prob_set_fp(mask_liq, :) = calculateFP( ...
        GFSetSimValues.setPGD(mask_liq, :), ...
        fra_media_set, fra_disper_set);
    prob_set_fp = prob_set_fp .* GFSetSimValues.liquProb;

    % -------------------- Landslide failure -----------------------------
    [N, M] = size(GFLandSimValues.landPGD);
    if isa(GFLandSimValues.landPGD, 'gpuArray')
        prob_land_fp = gpuArray.zeros(N, M, 'single');
    else
        prob_land_fp = zeros(N, M);
    end
    mask_land = sum(GFLandSimValues.landProb > 0.01, 2) > 0;
    prob_land_fp(mask_land, :) = calculateFP( ...
        GFLandSimValues.landPGD(mask_land, :), ...
        fra_media_land, fra_disper_land);
    prob_land_fp = prob_land_fp .* GFLandSimValues.landProb;

    % -------------------- Zone aggregation ------------------------------
    prob_lat_fp(prob_lat_fp   < 1e-6) = single(0);
    prob_set_fp(prob_set_fp   < 1e-6) = single(0);
    prob_land_fp(prob_land_fp < 1e-6) = single(0);

    [unique_zones, ~, bid_zid] = unique( ...
        [GFLatSimValues.bid_sid(:, end-1:end), ...
         GFSetSimValues.bid_sid(:, end), ...
         GFLandSimValues.bid_sid(:, end)], ...
        'rows');

    buildingGFfpEx.zone_GFfp_ex = max( ...
        prob_lat_fp(unique_zones(:, 2), :), ...
        prob_set_fp(unique_zones(:, 3), :));

    buildingGFfpCp.zone_GFfp_cp = max( ...
        0.2 * buildingGFfpEx.zone_GFfp_ex, ...
        prob_land_fp(unique_zones(:, 4), :));

    buildingGFfpEx.bid_zid = bid_zid;
    buildingGFfpCp.bid_zid = bid_zid;
end
