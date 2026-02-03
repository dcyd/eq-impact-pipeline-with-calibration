function [regionLossPerScenario, buildingLossRateScenarioSum] = ...
    calculate_build_structural_loss(sim_dir, stru_cost, batchSize)
%CALCULATE_BUILD_STRUCTURAL_LOSS Batched structural loss accumulation.
%   [regionLossPerScenario, buildingLossRateScenarioSum] =
%       CALCULATE_BUILD_STRUCTURAL_LOSS(sim_dir, stru_cost, batchSize)
%   reads pre-computed building damage probabilities for four damage
%   states (slight, moderate, extensive, complete) and aggregates expected
%   structural losses per scenario and per building in a batched manner.
%
%   The function:
%       • Phase 1 (COMP): computes complete-damage loss directly.
%       • Phase 2 (EXTE→MODE→SLIG): computes incremental losses for
%         extensive, moderate, and slight damage by differencing with the
%         previous higher state.
%
%   INPUTS
%       sim_dir   - Directory where *.mat files for each damage state are
%                   stored. Each file must contain:
%                       comp.mat, exte.mat, mode.mat, slig.mat
%                   with variables:
%                       buildingGSfp   - struct with fields:
%                           .unique_GSfp (N_unique×N_sce)   % GS probs
%                           .gs_bid      (N_B×1)            % mapping
%                       buildingGFfp   - either:
%                           (a) struct with fields
%                               .zone_GFfp_cp or .zone_GFfp_ex (N_zone×N_sce)
%                               .bid_zid (N_B×1)   % mapping building→zone
%                           (b) numeric matrix (N_B×N_sce) of GF probs
%
%       stru_cost - N_B×1 vector of structural replacement cost per
%                   building (same ordering as in gs_bid / bid_zid).
%
%       batchSize - (optional) number of scenarios processed per batch.
%                   Default: 50.
%
%   OUTPUTS
%       regionLossPerScenario - struct with fields:
%           .lossGS   - 1×N_sce expected loss from ground shaking only.
%           .lossGSGF - 1×N_sce expected loss from combined shaking and
%                       ground failure.
%
%       buildingLossRateScenarioSum - struct with fields:
%           .lossRateGS    - N_B×1 total (over scenarios) loss rate due to
%                            ground shaking only.
%           .lossRateGSGF  - N_B×1 total (over scenarios) loss rate due to
%                            combined shaking and ground failure.
%           .N_scenarios   - scalar, number of scenarios N_sce.
%
%   EXAMPLE
%       sim_dir   = 'results/NPT_demo/';
%       stru_cost = building_cost_vector;  % N_B×1
%       [regionLoss, bldLoss] = calculate_build_structural_loss( ...
%           sim_dir, stru_cost, 100);
%
%   SEE ALSO
%       calculate_building_damage_prob,
%       generate_building_quick_sim_data

    % ------------------------ Defaults -----------------------------------
    if nargin < 3 || isempty(batchSize)
        batchSize = 50;
    end

    damageStates = {'slig','mode','exte','comp'};
    % Expected loss ratios for [slight, moderate, extensive, complete]
    ELR = single([0.05, 0.20, 0.60, 1.00]);

    % =====================================================================
    % Phase 1: COMPLETE DAMAGE (baseline)
    % =====================================================================
    load(strcat(sim_dir, 'comp.mat'), 'buildingGSfp', 'buildingGFfp');

    uniBuildingGSfp = buildingGSfp.unique_GSfp;   % N_unique×N_sce
    gsBid           = buildingGSfp.gs_bid;        % N_B×1
    clear buildingGSfp

    % Size bookkeeping
    N_B   = size(gsBid, 1);          % number of buildings
    N_sce = size(uniBuildingGSfp, 2);% number of scenarios
    nbatch = ceil(N_sce / batchSize);

    % Result buffers
    region_GS_loss    = zeros(1, N_sce, 'single');
    region_GSGF_loss  = zeros(1, N_sce, 'single');
    build_GS_LR_sum   = zeros(N_B, 1, 'single');
    build_GSGF_LR_sum = zeros(N_B, 1, 'single');

    % Ground-failure representation
    isGFstruct = isstruct(buildingGFfp);
    if ~isGFstruct
        % already scenario-by-building matrix
        buildingGFfp = single(buildingGFfp);
    end

    % -------------------- Batch loop: complete state ---------------------
    for b = 1:nbatch
        idx = (b - 1) * batchSize + 1 : min(b * batchSize, N_sce);

        GSb = uniBuildingGSfp(gsBid, idx); % N_B×|idx|

        if isGFstruct
            GFb = single( ...
                buildingGFfp.zone_GFfp_cp(buildingGFfp.bid_zid, idx));
        else
            GFb = buildingGFfp(:, idx);
        end

        % Complete-state loss rate (LR) from GS and GS+GF
        LR_GS    = GSb * ELR(4);
        mix_comp = GSb + GFb - GSb .* GFb;
        LR_GSGF  = ELR(4) * mix_comp;

        % Building-level LR sums across scenarios
        build_GS_LR_sum   = build_GS_LR_sum   + sum(LR_GS,   2, 'native');
        build_GSGF_LR_sum = build_GSGF_LR_sum + sum(LR_GSGF, 2, 'native');

        % Convert LR to monetary loss
        LR_GS   = LR_GS   .* stru_cost;
        LR_GSGF = LR_GSGF .* stru_cost;

        % Regional aggregation per scenario
        region_GS_loss(idx)   = sum(LR_GS,   1, 'native');
        region_GSGF_loss(idx) = sum(LR_GSGF, 1, 'native');
    end

    % Store previous-level probabilities for incremental stages
    GS_prev_all    = uniBuildingGSfp;    % full unique GS matrix
    gs_bid_prev    = gsBid;             % mapping for previous level
    GF_prev_struct = buildingGFfp;      % may be struct or matrix

    clear GSb GFb LR_GS mix_comp LR_GSGF

    % =====================================================================
    % Phase 2: EXTE → MODE → SLIG (incremental losses)
    % =====================================================================
    for dsi = 3:-1:1   % 3=exte, 2=mode, 1=slig
        % Load current damage state probabilities
        matFile = strcat(sim_dir, [damageStates{dsi} '.mat']);
        load(matFile, 'buildingGSfp', 'buildingGFfp');

        uniBuildingGSfp = buildingGSfp.unique_GSfp;
        gsBid           = buildingGSfp.gs_bid;
        clear buildingGSfp

        GF_curr_struct = buildingGFfp;
        clear buildingGFfp

        if ~isGFstruct
            % If GF was a matrix in comp.mat, treat all levels as matrices
            GF_curr_all = single(GF_curr_struct);
        end

        for b = 1:nbatch
            idx = (b - 1) * batchSize + 1 : min(b * batchSize, N_sce);

            % ---- Current & previous ground-shaking slices ----
            GS_new = uniBuildingGSfp(gsBid, idx);
            GS_old = GS_prev_all(gs_bid_prev, idx);

            % ---- Current ground-failure (GF) slice ----
            if ~isempty(GF_curr_struct)
                if isGFstruct
                    if isfield(GF_curr_struct, "zone_GFfp_cp")
                        GF_new = single(GF_curr_struct.zone_GFfp_cp( ...
                            GF_curr_struct.bid_zid, idx));
                    elseif isfield(GF_curr_struct, "zone_GFfp_ex")
                        GF_new = single(GF_curr_struct.zone_GFfp_ex( ...
                            GF_curr_struct.bid_zid, idx));
                    else
                        error('GF_curr_struct must have zone_GFfp_cp or zone_GFfp_ex.');
                    end
                else
                    GF_new = GF_curr_all(:, idx);
                end

                mix_new = GS_new + GF_new - GS_new .* GF_new;
                clear GF_new
            else
                mix_new = GS_new;
            end

            % ---- Previous ground-failure (GF) slice ----
            if ~isempty(GF_prev_struct)
                if isGFstruct
                    if isfield(GF_prev_struct, "zone_GFfp_cp")
                        GF_old = single(GF_prev_struct.zone_GFfp_cp( ...
                            GF_prev_struct.bid_zid, idx));
                    elseif isfield(GF_prev_struct, "zone_GFfp_ex")
                        GF_old = single(GF_prev_struct.zone_GFfp_ex( ...
                            GF_prev_struct.bid_zid, idx));
                    else
                        error('GF_prev_struct must have zone_GFfp_cp or zone_GFfp_ex.');
                    end
                else
                    GF_old = GF_prev_struct(:, idx);
                end

                mix_old = GS_old + GF_old - GS_old .* GF_old;
                clear GF_old
            else
                mix_old = GS_old;
            end

            % ---- Incremental loss rates for this damage state ----
            LR_GS   = (GS_new - GS_old) * ELR(dsi);
            LR_GSGF = ELR(dsi) * (mix_new - mix_old);
            clear GS_new GS_old mix_new mix_old

            % Accumulate building-level LR sums
            build_GS_LR_sum   = build_GS_LR_sum   + sum(LR_GS,   2, 'native');
            build_GSGF_LR_sum = build_GSGF_LR_sum + sum(LR_GSGF, 2, 'native');

            % Convert LR to monetary loss
            LR_GS   = LR_GS   .* stru_cost;
            LR_GSGF = LR_GSGF .* stru_cost;

            % Regional aggregation per scenario
            region_GS_loss(idx)   = region_GS_loss(idx)   + sum(LR_GS,   1, 'native');
            region_GSGF_loss(idx) = region_GSGF_loss(idx) + sum(LR_GSGF, 1, 'native');

            clear LR_GS LR_GSGF
        end

        % ---- Rollover: current level becomes "previous" for next loop ----
        GS_prev_all    = uniBuildingGSfp;
        gs_bid_prev    = gsBid;
        GF_prev_struct = GF_curr_struct;

        clear uniBuildingGSfp GF_curr_struct
    end

    % =====================================================================
    % Pack outputs
    % =====================================================================
    regionLossPerScenario.lossGS   = region_GS_loss;
    regionLossPerScenario.lossGSGF = region_GSGF_loss;

    buildingLossRateScenarioSum.lossRateGS    = build_GS_LR_sum;
    buildingLossRateScenarioSum.lossRateGSGF  = build_GSGF_LR_sum;
    buildingLossRateScenarioSum.N_scenarios   = N_sce;
end
