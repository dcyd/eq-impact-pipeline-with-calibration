function [regionLossPerScenario, buildingLossRateScenarioSum] = ...
    calculate_build_other_loss(GSSimValues, building_vulnerability, cost_nonstr_cont, batch_size)
%CALCULATE_BUILD_OTHER_LOSS Non-structural and contents loss aggregation.
%   [regionLossPerScenario, buildingLossRateScenarioSum] =
%       CALCULATE_BUILD_OTHER_LOSS(GSSimValues, building_vulnerability, ...
%                                  cost_nonstr_cont, batch_size)
%   computes expected non-structural and contents loss per scenario and per
%   building, using vulnerability curves defined on intensity measure
%   levels (IMLs). The computation is batched over scenarios and uses the
%   GPU when available.
%
%   INPUTS
%       GSSimValues          - Structure with ground-shaking results for
%                              each building:
%                                .gmf_value : N_B×N_sce matrix of IM
%                                              (e.g., PGA or Sa) values
%                                .bid_sid   : N_B×K matrix; the last column
%                                              is a 1..N_B index used here.
%
%       building_vulnerability - Structure with vulnerability information:
%                                .vul_nonstrc_meanLRs : N_B×N_IML matrix of
%                                    mean loss ratios for non-structural
%                                    components at each IML.
%                                .vul_con_meanLRs     : N_B×N_IML matrix of
%                                    mean loss ratios for contents at
%                                    each IML.
%                                .vul_imls            : 1×N_IML or
%                                    N_IML×1 vector of IML values (must be
%                                    strictly increasing).
%
%       cost_nonstr_cont     - N_B×2 matrix of replacement costs:
%                                cost_nonstr_cont(:,1) = non-structural cost
%                                cost_nonstr_cont(:,2) = contents cost
%
%       batch_size           - (optional) batch size for processing
%                              scenarios. Default: 50.
%
%   OUTPUTS
%       regionLossPerScenario - struct with fields:
%           .lossNonStr   - 1×N_sce expected non-structural loss per
%                           scenario.
%           .lossCont     - 1×N_sce expected contents loss per scenario.
%
%       buildingLossRateScenarioSum - struct with fields:
%           .lossNonStr   - N_B×1 total (over scenarios) non-structural
%                           loss rate per building.
%           .lossCont     - N_B×1 total (over scenarios) contents loss
%                           rate per building.
%           .N_scenarios  - scalar, number of scenarios N_sce.
%
%   EXAMPLE
%       [regionLoss, bldLoss] = calculate_build_other_loss( ...
%           GSSimValues, building_vulnerability, cost_nonstr_cont, 100);
%
%   SEE ALSO
%       calculate_build_structural_loss, calculate_building_damage_prob

    % ------------------------ Defaults -----------------------------------
    if nargin < 4 || isempty(batch_size)
        batch_size = 50;
    end

    % ------------------------ Vulnerability inputs -----------------------
    vul_nonstrc_meanLRs = single(building_vulnerability.vul_nonstrc_meanLRs);
    vul_con_meanLRs     = single(building_vulnerability.vul_con_meanLRs);
    vul_imls            = single(building_vulnerability.vul_imls(1, :));
    clear building_vulnerability

    assert(issorted(vul_imls, 'strictascend'), ...
        'vul_imls must be strictly increasing.');

    % IML bin edges (midpoints between IMLs, with -inf/+inf at ends)
    edges = [-inf, ...
             (vul_imls(1:end-1) + vul_imls(2:end)) / 2, ...
              +inf];

    % ------------------------ Ground-motion inputs ------------------------
    bid_idx      = GSSimValues.bid_sid(:, end);    % N_B×1
    gmf_bid_zone = single(GSSimValues.gmf_value);  % N_B×N_sce (CPU)
    clear GSSimValues

    N_B   = size(bid_idx, 1);
    N_sce = size(gmf_bid_zone, 2);
    nbatch = ceil(N_sce / batch_size);

    % ------------------------ GPU / CPU setup ----------------------------
    try
        g      = gpuDevice();         % If no GPU, this will throw
        useGPU = g.DeviceSupported;

        rowIdx_GPU = single(gpuArray((1:N_B).'));

        vul_nonstrc_meanLRs = gpuArray(vul_nonstrc_meanLRs);
        vul_con_meanLRs     = gpuArray(vul_con_meanLRs);
        edges               = gpuArray(edges);
        gmf_bid_zone        = gpuArray(gmf_bid_zone);

    catch
        useGPU   = false;
        rowIdx_GPU = (1:N_B).';
    end

    % ------------------------ Result buffers -----------------------------
    region_nonStr_loss = zeros(1, N_sce, 'single');
    region_cont_loss   = zeros(1, N_sce, 'single');

    build_nonStr_LR_sum = zeros(N_B, 1, 'single');
    build_cont_LR_sum   = zeros(N_B, 1, 'single');

    % =====================================================================
    % Batched loss calculation
    % =====================================================================
    for b = 1:nbatch
        idx    = (b - 1) * batch_size + 1 : min(b * batch_size, N_sce);
        % i_nsce = numel(idx);  % not used but can be kept for debugging

        if useGPU
            % -------------------- GPU branch -----------------------------
            % Nb×|idx| bin indices over edges
            idx_gpu = discretize(gmf_bid_zone(bid_idx, idx), edges);

            % Linear indices into N_B×N_IML tables: (col-1)*N_B + row
            linIdx_gpu = (idx_gpu - 1) .* N_B + rowIdx_GPU;
            clear idx_gpu

            % ---- Non-structural loss rate ----
            LR_non_str_cols = vul_nonstrc_meanLRs(linIdx_gpu);
            build_nonStr_LR_sum = build_nonStr_LR_sum + ...
                                  sum(LR_non_str_cols, 2, 'native');

            LR_non_str_cols = LR_non_str_cols .* cost_nonstr_cont(:, 1);
            region_nonStr_loss(idx) = sum(LR_non_str_cols, 1, 'native');
            clear LR_non_str_cols

            % ---- Contents loss rate ----
            LR_con_cols = vul_con_meanLRs(linIdx_gpu);
            clear linIdx_gpu

            build_cont_LR_sum = build_cont_LR_sum + ...
                                sum(LR_con_cols, 2, 'native');

            LR_con_cols = LR_con_cols .* cost_nonstr_cont(:, 2);
            region_cont_loss(idx) = sum(LR_con_cols, 1, 'native');
            clear LR_con_cols

        else
            % -------------------- CPU branch -----------------------------
            % Nb×|idx| bin indices over edges
            idx_cpu = discretize(gmf_bid_zone(bid_idx, idx), edges);

            % Linear indices into N_B×N_IML tables: (col-1)*N_B + row
            linIdx_cpu = (idx_cpu - 1) .* N_B + rowIdx_GPU;
            clear idx_cpu

            % ---- Non-structural loss rate ----
            LR_non_str_cols = vul_nonstrc_meanLRs(linIdx_cpu);
            build_nonStr_LR_sum = build_nonStr_LR_sum + ...
                                  sum(LR_non_str_cols, 2, 'native');

            LR_non_str_cols = LR_non_str_cols .* cost_nonstr_cont(:, 1);
            region_nonStr_loss(idx) = sum(LR_non_str_cols, 1, 'native');
            clear LR_non_str_cols

            % ---- Contents loss rate ----
            LR_con_cols = vul_con_meanLRs(linIdx_cpu);
            clear linIdx_cpu

            build_cont_LR_sum = build_cont_LR_sum + ...
                                sum(LR_con_cols, 2, 'native');

            LR_con_cols = LR_con_cols .* cost_nonstr_cont(:, 2);
            region_cont_loss(idx) = sum(LR_con_cols, 1, 'native');
            clear LR_con_cols
        end
    end

    % =====================================================================
    % Pack outputs
    % =====================================================================
    regionLossPerScenario.lossNonStr = region_nonStr_loss;
    regionLossPerScenario.lossCont   = region_cont_loss;

    buildingLossRateScenarioSum.lossNonStr  = build_nonStr_LR_sum;
    buildingLossRateScenarioSum.lossCont    = build_cont_LR_sum;
    buildingLossRateScenarioSum.N_scenarios = N_sce;
end
