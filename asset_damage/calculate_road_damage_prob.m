function [bridge_GSfp,edge_Block_GSfp,edge_GFfp,edge_Block_GFfp] = calculate_road_damage_prob(BuildingInfo, RoadSystem, FreeParams)

    % GSSimValues, GFLatSimValues, GFSetSimValues,GFLandSimValues,
% SeismicScenario.PGA: % N_Z x (3+N_sce) matrix, where N_Z is the number of zones, the first three columns are the zone ID, longitude,  and latitude, and the remaining columns are the PGA values for different scenarios.
% SeismicScenario.SA0p3: same as PGA, but for SA(0.3s)
% SeismicScenario.SA0p6: same as PGA, but for SA(0.6s) 
% SeismicScenario.SA1p0: same as PGA, but for SA(1.0s)
% SeismicScenario.landPGD: same as PGA, but for land PGD
% SeismicScenario.latPGD: same as PGA, but for lateral PGD
% SeismicScenario.setPGD: same as PGA, but for settlement PGD

% BuildingInfo.building_data: % N_B x 3 matrix, where N_B is the number of buildings, the first three columns are the building ID, longitude, latitude
% BuildingInfo.building_zone: % N_B x 2 matrix, the building ID and the zone ID for each building.
% BuildingInfo.building_edge_data: % N_B x 6 matrix, where the first column is the building ID, the second column is the edge ID, and the remaining columns are other attributes.
% BuildingInfo.fra_params: % A structure containing the parameters for the fragility curves, including:
%  .fra_inten_code: % N_B x 1 matrix, the intensity measure code for each building. 1: PGA, 2: SA(0.3s), 3: SA(0.6s), 4: SA(1.0s).
%  .fra_media: % N_B x Ns matrix, the median values of fragility curves for Ns states. If size is 1 x Ns, the same fragility curves are applied to all points.
%  .fra_disper: % N_B x 1, 1 x Ns, N_B x 1, or N_B x Ns matrix, the dispersion values of fragility curves. If 1 x 1, the same dispersion is applied to all states and points. If 1 x Ns, the same dispersion is applied to all points but varies by state. If N_B x 1, the same dispersion is applied to all states but varies by point. If N_B x Ns, dispersion varies by both state and point.
% BuildingInfo.wd_mean_std: % N_B_E x 2 matrix, where the first column is the mean value of the width and the second column is the standard deviation.

% RoadSystem.edge_data: % N_E x 8 matrix, where N_E is the number of edges, the first column is the edge ID, and the remaining columns are other attributes.
% RoadSystem.edge_zone: % N_E_Z x 2 matrix, where the first column is the edge ID and the second column is the zone ID.

% FreeParams: % A structure containing the parameters for the fragility curves and other calculations.
%  .cfBridgeGSMedia: Correction Factor of Median for bridge fragility function in ground shaking
% .cfUrbanRoadGFMedia: Correction Factor of Median for urban road fragility function in ground failure
% .cfMajorRoadGFMedia: Correction Factor of Median for major road fragility function in ground failure
% .cfBridgeGFMedia: Correction Factor of Median for bridge fragility function in ground failure
% .cfTunnelGFMedia: Correction Factor of Median for tunnel fragility function in ground failure
% .cfWDMean: Correction Factor of Mean for building debris width function
% .cfWDSTD: Correction Factor of Standard Deviation for building debris width function

% Params: % A structure containing additional parameters for the simulation.
% .N_sim: Number of simulations to perform (default is 1).
% .ses_range: Range of seismic scenarios to consider (default is all scenarios).

    %%
    if nargin < 3
        FreeParams = struct;
    end
    if isnumeric(FreeParams)
        buff = FreeParams;
        FreeParams = struct;
        FreeParams.cfBuildingLatMedia = buff(1);
        FreeParams.cfBuildingSetMedia = buff(2);
        FreeParams.cfBuildingLandMedia = buff(3);
        FreeParams.cfBuildGSMedia = buff(4);
        FreeParams.cfBridgeGSMedia = buff(5);
        FreeParams.cfUrbanRoadGFMedia = buff(6);
        FreeParams.cfMajorRoadGFMedia = buff(7);
        FreeParams.cfBridgeGFMedia = buff(8);
        FreeParams.cfTunnelGFMedia = buff(9);
        FreeParams.cfWDMean = buff(10);
        FreeParams.cfWDSTD = buff(11);
    end


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
    if ~isfield(FreeParams, 'cfBridgeGSMedia')
        FreeParams.cfBridgeGSMedia = 1;
    end
    if ~isfield(FreeParams, 'cfUrbanRoadGFMedia')
        FreeParams.cfUrbanRoadGFMedia = 1;
    end
    if ~isfield(FreeParams, 'cfMajorRoadGFMedia')
        FreeParams.cfMajorRoadGFMedia = 1;
    end
    if ~isfield(FreeParams, 'cfBridgeGFMedia')
        FreeParams.cfBridgeGFMedia = 1;
    end
    if ~isfield(FreeParams, 'cfTunnelGFMedia')
        FreeParams.cfTunnelGFMedia = 1;
    end
    if ~isfield(FreeParams, 'cfWDMean')
        FreeParams.cfWDMean = 1;
    end
    if ~isfield(FreeParams, 'cfWDSTD')
        FreeParams.cfWDSTD = 1;
    end

    %%
    N_sce = size(RoadSystem.GSRoadSimValues.SA1p0,2);
    N_E = size(RoadSystem.GSRoadSimValues.is_bri,1);
    isgpu = isa(RoadSystem.GSRoadSimValues.SA1p0, 'gpuArray');

    %% GS
    if isgpu
        bridge_GSfp = gpuArray.zeros(N_E,N_sce,'single');
    else
        bridge_GSfp = zeros(N_E,N_sce);
    end
    is_bri = RoadSystem.GSRoadSimValues.is_bri;
    if any(is_bri)
        bridge_GSfp(is_bri,:) = calculateBridgeGSfp(RoadSystem.GSRoadSimValues, FreeParams);  % P3
    end
    %% GF

    edge_GFfp = calculateRoadGFfp(RoadSystem.LandRoadSimValues, RoadSystem.LiqRoadSimValues, FreeParams);  % P4

    %% building debris
    [edge_Block_GSfp,edge_Block_GFfp] = calculateRoadBlockFp(BuildingInfo, RoadSystem.edge_id_type(:,1),FreeParams);  % P5
    
    %% unique
    try
        bridge_GSfp(bridge_GSfp<10^(-6))=single(0);
        [uni_bridge_GSfp,~,gs_bid] = unique(bridge_GSfp,"rows");
    catch
        bridge_GSfp = gather(bridge_GSfp);
        [uni_bridge_GSfp,~,gs_bid] = unique(bridge_GSfp,"rows");
    end
    clear bridge_GSfp
    bridge_GSfp.unique_GSfp = uni_bridge_GSfp;
    bridge_GSfp.gs_bid = gs_bid;
    
    try
        edge_GFfp(edge_GFfp<10^(-6))=single(0);
        [uni_edge_GFfp,~,gf_eid] = unique(edge_GFfp,"rows");
    catch
        edge_GFfp = gather(edge_GFfp);
        [uni_edge_GFfp,~,gf_eid] = unique(edge_GFfp,"rows");
    end
    clear edge_GFfp
    edge_GFfp.unique_GFfp = uni_edge_GFfp;
    edge_GFfp.gf_eid = gf_eid;
    
    try
        edge_Block_GSfp(edge_Block_GSfp<10^(-6))=single(0);
        [uni_edge_blockGSfp,~,gs_block_eid] = unique(edge_Block_GSfp,"rows");
    catch
        edge_Block_GSfp = gather(edge_Block_GSfp);
        [uni_edge_blockGSfp,~,gs_block_eid] = unique(edge_Block_GSfp,"rows");
    end
    clear edge_Block_GSfp
    edge_Block_GSfp.unique_blockGSfp = uni_edge_blockGSfp;
    edge_Block_GSfp.gs_block_eid = gs_block_eid;
    
    try
        edge_Block_GFfp(edge_Block_GFfp<10^(-6))=single(0);
        [uni_edge_blockGFfp,~,gf_block_eid] = unique(edge_Block_GFfp,"rows");
    catch
        edge_Block_GFfp = gather(edge_Block_GFfp);
        [uni_edge_blockGFfp,~,gf_block_eid] = unique(edge_Block_GFfp,"rows");
    end
    clear edge_Block_GFfp
    edge_Block_GFfp.unique_blockGFfp = uni_edge_blockGFfp;
    edge_Block_GFfp.gf_block_eid = gf_block_eid;

    %% Combine
        % building_GFfp = calculateBuildGFfp(build_gmf,FreeParams);  % P2


    % 
    % buildingFP = building_GFfp+building_GSfp - building_GFfp.*building_GSfp;
    % roadFP = 1-(1-bridge_GSfp).*(1-edge_Block_GSfp).*(1-max(edge_GFfp,edge_Block_GFfp));
end




function prob = calculateFP(gmf_value, fra_media, fra_disper)
    is_batch = 0;
    batch_size = 800; % 可根据内存情况调整

    % 批处理参数
    [N, N_sce] = size(gmf_value);
    Ns = size(fra_media, 2);

    % 扩展参数
    if size(fra_media, 1) == 1
        fra_media = repmat(fra_media, N, 1);
    end
    if isscalar(fra_disper)
        fra_disper = repmat(fra_disper, N, Ns);
    elseif size(fra_disper, 1) == 1
        fra_disper = repmat(fra_disper, N, 1);
    elseif size(fra_disper, 2) == 1
        fra_disper = repmat(fra_disper, 1, Ns);
    end

    isgpu = isa(gmf_value, 'gpuArray');
    if isgpu
        if ~isa(fra_media, 'gpuArray')
            fra_media = single(gpuArray(fra_media));
            fra_disper = single(gpuArray(fra_disper));
        end
    end


    fra_media_exp  = reshape(fra_media,  N, 1, Ns);   % N×1×Ns
    fra_disper_exp = reshape(fra_disper, N, 1, Ns);   % N×1×Ns
    clear fra_media fra_disper
    log_fra_media_exp = log(fra_media_exp);
    clear fra_media_exp
    threshold   = log_fra_media_exp - 3*fra_disper_exp;
    mask = sum(log(gmf_value) > threshold,2)>0;
    clear threshold
    if isgpu
        prob = gpuArray.zeros(N, N_sce, Ns, 'single');
    else
        prob = zeros(N, N_sce, Ns);
    end
    if is_batch
        num_batches = ceil(N_sce / batch_size);
        for b = 1:num_batches
            idx_start = (b-1)*batch_size + 1;
            idx_end   = min(b*batch_size, N_sce);
            idx       = idx_start:idx_end;
            % prob(mask, idx, :) = normcdf( numerator(mask,idx) ./ (fra_disper_exp(mask,:) + eps) );
            prob(mask, idx, :) = logncdf(gmf_value(mask,idx),log_fra_media_exp(mask,:),fra_disper_exp(mask,:));
        end
    else
        prob(mask, :, :) = logncdf(gmf_value(mask,:),log_fra_media_exp(mask,:),fra_disper_exp(mask,:));
    end
    clear fra_disper_exp threshold mask gmf_value log_fra_media_exp

end




function bridge_GSfp = calculateBridgeGSfp(GSRoadSimValues, params)

    if nargin<2
        params = struct(); % Initialize params if not provided
    end
    % Extract or set default parameters
    if isfield(params, 'bri_gs_fra_media')
        fra_media = params.bri_gs_fra_media;  % Median for settlement
    else
        fra_media = 1.1;  % Default median as 1.1 from Hazus (the middle one)
        % fprintf('Default median value %f is used for bridge fragility curves.\n', fra_media);
    end
    if isfield(params, 'bri_gs_fra_disper')
        fra_disper = params.bri_gs_fra_disper;  % Dispersion for settlement
    else
        fra_disper = 0.6; % Default dispersion as 0.6 from Hazus
        % fprintf('Default dispersion value %f is used for bridge fragility curves.\n', fra_disper);
    end

    if isfield(params, 'cfBridgeGSMedia')
        cfBridgeGSMedia = params.cfBridgeGSMedia;
    else
        cfBridgeGSMedia = 1;
    end

    fra_media = fra_media*cfBridgeGSMedia;  % Adjust median values for building height

    zi_bridge_GSfp = calculateFP(GSRoadSimValues.SA1p0, fra_media, fra_disper);
    bridge_GSfp = zi_bridge_GSfp(GSRoadSimValues.bri_sids(:,2),:);
end

function tunnel_GSfp = calculateTunnelGSfp(tunnel_gmf, params)

% tunnel_gmf An N_B, N_sce, containing PGA
    nd_bri = ndims(tunnel_gmf);
    if nd_bri < 3
        [N, M] = size(tunnel_gmf);
        tunnel_gmf = reshape(tunnel_gmf, [N, 1, M]);
    end
    if nargin<2
        params = struct(); % Initialize params if not provided
    end
    % Extract or set default parameters
    if isfield(params, 'fra_media')
        fra_media = params.fra_media;  % Median for settlement
    else
        fra_media = 1.1;  % Default median as 1.1 from Hazus (the middle one)
        % fprintf('Default median value %f is used for tunnel fragility curves.\n', fra_media);
    end
    if isfield(params, 'fra_disper')
        fra_disper = params.fra_disper;  % Dispersion for settlement
    else
        fra_disper = 0.6; % Default dispersion as 0.6 from Hazus
        % fprintf('Default dispersion value %f is used for tunnel fragility curves.\n', fra_disper);
    end

    if isfield(params, 'cf_tun_media')
        k_bf = params.cf_tun_media;
    else
        k_bf = 1;
    end

    fra_media = fra_media*k_bf;  % Adjust median values for building height

    tunnel_GSfp = calculateFP(tunnel_gmf, fra_media, fra_disper);
end


function [road_GFfp] = calculateRoadGFfp(LandRoadSimValues, LiqRoadSimValues, params)

    isgpu = isa(LandRoadSimValues.landPGD, 'gpuArray');
    if isfield(params, 'cfBridgeGFMedia')
        cfBridgeGFMedia = params.cfBridgeGFMedia;
    else
        cfBridgeGFMedia = 1;
    end
    if isfield(params, 'cfTunnelGFMedia')
        cfTunnelGFMedia = params.cfTunnelGFMedia;
    else
        cfTunnelGFMedia = 1;
    end
    if isfield(params, 'cfUrbanRoadGFMedia')
        cfUrbanRoadGFMedia = params.cfUrbanRoadGFMedia;
    else
        cfUrbanRoadGFMedia = 1;
    end
    if isfield(params, 'cfUrbanRoadGFMedia')
        cfMajorRoadGFMedia = params.cfMajorRoadGFMedia;
    else
        cfMajorRoadGFMedia = 1;
    end

    function [fra_media,fra_disper] = sample_road_gf_fra(edge_type)
        % The fragility functions from HAZUS Manuel
        bri_media = [13.8]*cfBridgeGFMedia;
        bri_disper = 0.2;
    
        tun_media = [60]*cfTunnelGFMedia;
        tun_disper = 0.5;
    
        maj_road_media = [60]*cfMajorRoadGFMedia;
        urb_road_media = [24]*cfUrbanRoadGFMedia;
        road_disper = 0.7;
    
        % Median and dispersion for each edge type
        % 1: major road; 2: urban road; 3: bridge; 4: tunnel
        is_e_maj = edge_type == 1;
        is_e_urb = edge_type == 2;
        is_bridge = edge_type == 3;
        is_tunnel = edge_type == 4;
        N_E = size(edge_type,1);
        if isgpu
            fra_media = gpuArray.zeros(N_E,1,'single');
            fra_disper = gpuArray.zeros(N_E,1,'single');
        else
            fra_media = zeros(N_E,1);
            fra_disper = zeros(N_E,1);
        end
        fra_media(is_e_maj,:) = repmat(maj_road_media, sum(is_e_maj), 1);
        fra_media(is_e_urb,:) = repmat(urb_road_media, sum(is_e_urb), 1);
        fra_disper(:) = road_disper;
    
        fra_media(is_bridge,:) = repmat(bri_media, sum(is_bridge), 1);
        fra_disper(is_bridge) = bri_disper;
        fra_media(is_tunnel,:) = repmat(tun_media, sum(is_tunnel), 1);
        fra_disper(is_tunnel) = tun_disper;
    end
    [fra_media,fra_disper] = sample_road_gf_fra(LiqRoadSimValues.edge_type);
    road_liq_GFfp = calculateFP(LiqRoadSimValues.liquPGD, fra_media, fra_disper);
    road_liq_GFfp = road_liq_GFfp.*LiqRoadSimValues.liquProb; %  max_liq_pro_per_edge
    [fra_media,fra_disper] = sample_road_gf_fra(LandRoadSimValues.edge_type);
    road_land_GFfp = calculateFP(LandRoadSimValues.landPGD, fra_media, fra_disper);
    road_land_GFfp = road_land_GFfp.*LandRoadSimValues.landProb; % max_land_pro_per_edge
    
    road_GFfp = max(road_liq_GFfp(LiqRoadSimValues.edge_id_type_sids(:,end),:),road_land_GFfp(LandRoadSimValues.edge_id_type_sids(:,end),:));
end



function [edge_Block_GSfp,edge_Block_GFfp] = calculateRoadBlockFp(BuildingInfo, edge_ids, params)

    if isfield(params, 'cfWDMean')
        cfWDMean = params.cfWDMean;
    else
        cfWDMean = 1;
    end
    if isfield(params, 'cfWDSTD')
        cfWDSTD = params.cfWDSTD;
    else
        cfWDSTD = 1;
    end
    if isfield(params, 'thres_width')
        thres_width = params.thres_width;
    else
        thres_width = 3.5;
    end

    building_edge_data = BuildingInfo.building_edge_data;
    wd_mean_std = BuildingInfo.wd_mean_std;
    wd_mean_std(:,1) = wd_mean_std(:,1)*cfWDMean;
    wd_mean_std(:,2) = wd_mean_std(:,2)*cfWDSTD;

    N_sce = size(BuildingInfo.GSSimValues.gmf_value,2);
    % Retrieve necessary values from `building_edge_data`
    blockage_threshold = building_edge_data(:,3) + building_edge_data(:,6) - thres_width;

    edge_build_fp_wd = 1-normcdf(blockage_threshold, wd_mean_std(:,1), wd_mean_std(:,2));

    [building_GSfp,building_GFfp] = calculate_building_damage_prob(BuildingInfo, params);
    [~,sel_build_tag] = ismember(BuildingInfo.building_edge_data(:,1),BuildingInfo.GSSimValues.bid_sid(:,1));
    edge_build_GS_fp = edge_build_fp_wd.*building_GSfp.unique_GSfp(building_GSfp.gs_bid(sel_build_tag),:);
    clear building_GSfp
    if isstruct(building_GFfp)
        edge_build_GF_fp = edge_build_fp_wd.*building_GFfp.zone_GFfp_cp(building_GFfp.bid_zid(sel_build_tag,:),:,end);
    else
        edge_build_GF_fp = edge_build_fp_wd.*building_GFfp(sel_build_tag,:,end);
    end
    clear building_GFfp
    edge_build_GS_fp = edge_build_GS_fp.*BuildingInfo.GSSimValues.rate_cd2collap(BuildingInfo.GSSimValues.bid_sid(sel_build_tag,end),:);% 乘以complete 2 collapsed的比例
    edge_build_GF_fp = edge_build_GF_fp.*BuildingInfo.GSSimValues.rate_cd2collap(BuildingInfo.GSSimValues.bid_sid(sel_build_tag,end),:);% 乘以complete 2 collapsed的比例
    clear BuildingInfo
    [~,edge_tag] = ismember(building_edge_data(:, 2), edge_ids);

    % edge_tag_expanded = repmat(edge_tag, 1, N_sce);

    % 使用 accumarray 计算所有列的结果
    % edge_tag_2 = repelem((1:N_sce)', numel(edge_tag));
    try
        edge_tag_expanded = repmat(edge_tag, 1, N_sce);
        edge_tag_2 = repelem((1:N_sce)', numel(edge_tag));
        edge_Block_GSfp = 1 - accumarray([edge_tag_expanded(:), edge_tag_2(:)], ...
                                         1 - edge_build_GS_fp(:), ...
                                         [numel(edge_ids), N_sce], ...
                                         @prod);
        
        edge_Block_GFfp = 1 - accumarray([edge_tag_expanded(:), edge_tag_2(:)], ...
                                         1 - edge_build_GF_fp(:), ...
                                         [numel(edge_ids), N_sce], ...
                                         @prod);
        clear edge_tag_expanded edge_tag_2
    catch ME
        batch_size = 50; % 每批处理的场景数，可根据内存调整
        num_batches = ceil(N_sce / batch_size);

        edge_Block_GSfp = zeros(numel(edge_ids), N_sce, 'like', edge_build_GS_fp);
        edge_Block_GFfp = zeros(numel(edge_ids), N_sce, 'like', edge_build_GF_fp);

        for b = 1:num_batches
            idx_start = (b-1)*batch_size + 1;
            idx_end = min(b*batch_size, N_sce);
            idx = idx_start:idx_end;

            % 批量取出
            edge_build_GS_fp_batch = edge_build_GS_fp(:, idx);
            edge_build_GF_fp_batch = edge_build_GF_fp(:, idx);

            edge_tag_2_batch = repelem((1:length(idx))', numel(edge_tag));
            edge_tag_expanded_batch = repmat(edge_tag, 1, length(idx));

            % GS
            edge_Block_GSfp(:, idx) = 1 - accumarray([edge_tag_expanded_batch(:), edge_tag_2_batch(:)], ...
                                                    1 - edge_build_GS_fp_batch(:), ...
                                                    [numel(edge_ids), length(idx)], ...
                                                    @prod);

            % GF
            edge_Block_GFfp(:, idx) = 1 - accumarray([edge_tag_expanded_batch(:), edge_tag_2_batch(:)], ...
                                                    1 - edge_build_GF_fp_batch(:), ...
                                                    [numel(edge_ids), length(idx)], ...
                                                    @prod);
        end
        % % 将 GPU 数据转移到 CPU 进行 accumarray 操作
        % if isa(edge_build_GS_fp, 'gpuArray')
        %     edge_build_GS_fp = gather(edge_build_GS_fp);
        %     edge_build_GF_fp = gather(edge_build_GF_fp);
        %     edge_tag_expanded = gather(edge_tag_expanded);
        %     edge_tag_2 = gather(edge_tag_2);
        % end
        
        % edge_Block_GSfp = 1 - accumarray([edge_tag_expanded(:), edge_tag_2(:)], ...
        %                                  1 - edge_build_GS_fp(:), ...
        %                                  [numel(edge_ids), N_sce], ...
        %                                  @prod);
        
        % edge_Block_GFfp = 1 - accumarray([edge_tag_expanded(:), edge_tag_2(:)], ...
        %                                  1 - edge_build_GF_fp(:), ...
        %                                  [numel(edge_ids), N_sce], ...
        %                                  @prod);
    end  
    
    % 将未匹配的 edge_tag 设置为 0
    unmatched_edges = ~ismember((1:numel(edge_ids))', edge_tag);
    edge_Block_GSfp(unmatched_edges, :) = 0;
    edge_Block_GFfp(unmatched_edges, :) = 0;

end

