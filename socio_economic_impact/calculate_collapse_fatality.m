function [regionBuildColl, regionFatality, buildFatalitySum] = ...
    calculate_collapse_fatality(building_GSfp, building_GFfp, buildingInfo, batch_size)
%CALCULATE_COLLAPSE_FATALITY Expected collapses and fatalities per scenario.
%   [regionBuildColl, regionFatality, buildFatalitySum] =
%       CALCULATE_COLLAPSE_FATALITY(building_GSfp, building_GFfp, ...
%                                   buildingInfo, batch_size)
%   computes expected numbers of collapsed buildings and fatalities for
%   each seismic scenario by combining ground-shaking and ground-failure
%   complete-damage probabilities with collapse and lethality rates.
%
%   INPUTS
%       building_GSfp - Structure with building-level ground-shaking
%           complete-damage probabilities:
%               .unique_GSfp - N_unique×N_sce matrix of complete-damage
%                              exceedance probabilities from shaking.
%               .gs_bid      - N_B×1 index mapping each building to a row
%                              in unique_GSfp.
%
%       building_GFfp - Structure with building-level ground-failure
%           complete-damage probabilities:
%               .zone_GFfp_cp - N_zone×N_sce matrix of complete-damage
%                               exceedance probabilities from ground
%                               failure (e.g., liquefaction, landslides).
%               .bid_zid      - N_B×1 index mapping each building to a
%                               zone row in zone_GFfp_cp.
%
%       buildingInfo  - Structure containing building attributes:
%               .building_pop_lethality - N_B×K matrix; column 7 is the
%                                         number of occupants per building.
%               .fra_params             - structure with:
%                   .rate_cd2collap    - N_B×1 ratio of complete damage
%                                        to collapse for each building.
%                   .rate_collap2fatali- N_B×1 rate of fatalities per
%                                        collapsed occupant.
%
%       batch_size    - (optional) number of scenarios processed per batch.
%                       Default: 50.
%
%   OUTPUTS
%       regionBuildColl - 1×N_sce vector of expected number of collapsed
%                         buildings per scenario (sum over all buildings).
%
%       regionFatality  - 1×N_sce vector of expected fatalities per
%                         scenario (sum over all buildings).
%
%       buildFatalitySum - N_B×1 vector of total expected fatalities per
%                          building aggregated over all scenarios.
%
%   EXAMPLE
%       [regColl, regFat, bldFat] = calculate_collapse_fatality( ...
%           building_GSfp, building_GFfp, buildingInfo, 100);
%
%   SEE ALSO
%       calculate_building_damage_prob, calculate_build_structural_loss

    % ------------------------ Defaults -----------------------------------
    if nargin < 4 || isempty(batch_size)
        batch_size = 50;
    end

    % ------------------------ Unpack inputs ------------------------------
    uni_building_GSfp = building_GSfp.unique_GSfp;    % N_unique×N_sce
    gs_bid            = building_GSfp.gs_bid;         % N_B×1

    zone_GFfp_cp = building_GFfp.zone_GFfp_cp;        % N_zone×N_sce
    gf_bid_zid   = building_GFfp.bid_zid;             % N_B×1

    clear building_GSfp building_GFfp

    building_occupants   = buildingInfo.building_pop_lethality(:, 7); % N_B×1
    rate_collap2fatali   = buildingInfo.fra_params.rate_collap2fatali; % N_B×1
    rate_cd2collap       = buildingInfo.fra_params.rate_cd2collap;     % N_B×1

    clear buildingInfo

    % ------------------------ Size bookkeeping ---------------------------
    N_B   = size(gs_bid, 1);               % number of buildings
    N_sce = size(uni_building_GSfp, 2);    % number of scenarios
    nbatch = ceil(N_sce / batch_size);

    % Result buffers
    regionBuildColl  = zeros(1, N_sce, 'single');
    regionFatality   = zeros(1, N_sce, 'single');
    buildFatalitySum = zeros(N_B, 1, 'single');

    % =====================================================================
    % Batched calculation of collapses and fatalities
    % =====================================================================
    for b = 1:nbatch
        idx = (b - 1) * batch_size + 1 : min(b * batch_size, N_sce);

        % Shaking-related complete-damage prob (per building, per scenario)
        GSb = uni_building_GSfp(gs_bid, idx);            % N_B×|idx|
        % Ground-failure-related complete-damage prob
        GFb = single(zone_GFfp_cp(gf_bid_zid, idx));     % N_B×|idx|

        % Combined complete-damage probability (GS + GF - intersection)
        mix_comp_prob = GSb + GFb - GSb .* GFb;          % N_B×|idx|

        % Convert complete-damage probability to collapse probability
        mix_comp_prob_collapsed = gather(mix_comp_prob .* rate_cd2collap); 
        % N_B×|idx|

        % Expected fatalities = collapse_prob × occupants × fatality_rate
        mix_comp_fatality = gather( ...
            mix_comp_prob_collapsed .* building_occupants .* rate_collap2fatali);
        % N_B×|idx|

        % Regional totals per scenario
        regionBuildColl(1, idx) = sum(mix_comp_prob_collapsed, 1);
        regionFatality(1, idx)  = sum(mix_comp_fatality,       1);

        % Building totals across scenarios
        buildFatalitySum = buildFatalitySum + sum(mix_comp_fatality, 2);
    end
end
