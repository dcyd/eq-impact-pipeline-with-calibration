function [BuildingInfo_quick,RoadSystem_quick] = generate_road_quick_sim_data(BuildingInfo,RoadSystem,SeismicScenario,isGpu)

    if nargin <= 3
        isGpu = false;
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


    % The BuildingInfo nearby roads
    sel_build_idx = ismember(BuildingInfo.building_data(:,1),unique(BuildingInfo.building_edge_data(:,1)));
    sel_BuildingInfo = BuildingInfo;
    clear BuildingInfo;
    sel_BuildingInfo.building_data = sel_BuildingInfo.building_data(sel_build_idx,:);
    sel_BuildingInfo.building_zone = sel_BuildingInfo.building_zone(sel_build_idx,:);
    sel_BuildingInfo.fra_params.fra_disper = sel_BuildingInfo.fra_params.fra_disper(sel_build_idx,:);
    sel_BuildingInfo.fra_params.fra_media = sel_BuildingInfo.fra_params.fra_media(sel_build_idx,:);
    sel_BuildingInfo.fra_params.fra_inten_code = sel_BuildingInfo.fra_params.fra_inten_code(sel_build_idx,:);
    sel_BuildingInfo.fra_params.rate_cd2collap = sel_BuildingInfo.fra_params.rate_cd2collap(sel_build_idx,:);

    sel_BuildingInfo.fra_params.fra_media = sel_BuildingInfo.fra_params.fra_media(:,end);
    N_B = length(sel_BuildingInfo.building_data(:,1));
    N_sce = length(SeismicScenario.landPGD(1,:))-3;
    N_sta = size(sel_BuildingInfo.fra_params.fra_media, 2); % Number of damage states

    [~, gmf_zids] = ismember(sel_BuildingInfo.building_zone(:, 2), SeismicScenario.PGA(:, 1));

    idx = 1:N_sce;
    fra_inten_code = sel_BuildingInfo.fra_params.fra_inten_code;  % Intensity measure codes
    fra_media = sel_BuildingInfo.fra_params.fra_media;  % Median
    fra_disper = sel_BuildingInfo.fra_params.fra_disper;  % Dispersion
    rate_cd2collap = sel_BuildingInfo.fra_params.rate_cd2collap;

    % Record the needed gmf value for the selected buildings
    if isa(SeismicScenario.PGA, 'gpuArray')
        gmf_value = gpuArray.zeros(N_B, N_sce, 'single');
    else
        gmf_value = zeros(N_B, N_sce,'single');
    end

    PGA = SeismicScenario.PGA;
    SA0p3 = SeismicScenario.SA0p3;
    SA0p6 = SeismicScenario.SA0p6;
    SA1p0 = SeismicScenario.SA1p0;
    landPGD = SeismicScenario.landPGD;
    latPGD = SeismicScenario.latPGD;
    setPGD = SeismicScenario.setPGD;
    liquProb = SeismicScenario.liquProb;
    landProb = SeismicScenario.landProb;
    clear SeismicScenario

    sel_bui = fra_inten_code == 1;
    gmf_value(sel_bui,:) = PGA(gmf_zids(sel_bui),idx+3);
    sel_bui = fra_inten_code == 2;
    gmf_value(sel_bui,:) = SA0p3(gmf_zids(sel_bui),idx+3);
    sel_bui = fra_inten_code == 3;
    gmf_value(sel_bui,:) = SA0p6(gmf_zids(sel_bui),idx+3);
    sel_bui = fra_inten_code == 4;
    gmf_value(sel_bui,:) = SA1p0(gmf_zids(sel_bui),idx+3);

    gmf_fraMed_fraDis = single([gmf_value, fra_media, fra_disper,rate_cd2collap]);
    
    [uni_values, ~, build_sids] = unique(gmf_fraMed_fraDis, 'rows');
    GSSimValues.gmf_value = uni_values(:, 1:end-3);
    GSSimValues.fra_media = uni_values(:, end-2);
    GSSimValues.fra_disper = uni_values(:, end-1);
    GSSimValues.rate_cd2collap = uni_values(:, end);
    GSSimValues.bid_sid = [sel_BuildingInfo.building_data, build_sids];

    clear gmf_value fra_media fra_disper fra_inten_code gmf_fraMed_fraDis uni_values
    clear PGA SA0p3 SA0p6

    % Reorganize ground failure related data

    [uni_values, ~, build_sids] = unique(single([latPGD(gmf_zids, idx + 3),liquProb(gmf_zids, idx + 3)]), 'rows');
    GFLatSimValues.latPGD = uni_values(:, 1:N_sce);
    GFLatSimValues.liquProb = uni_values(:, N_sce+1:end);
    GFLatSimValues.bid_sid = [sel_BuildingInfo.building_data, build_sids];
    clear uni_values

    [uni_values, ~, build_sids] = unique(single([setPGD(gmf_zids, idx + 3),liquProb(gmf_zids, idx + 3)]), 'rows');
    GFSetSimValues.setPGD = uni_values(:, 1:N_sce);
    GFSetSimValues.liquProb = uni_values(:, N_sce+1:end);
    GFSetSimValues.bid_sid = [sel_BuildingInfo.building_data, build_sids];
    clear uni_values

    [uni_values, ~, build_sids] = unique(single([landPGD(gmf_zids, idx + 3),landProb(gmf_zids, idx + 3)]), 'rows');
    GFLandSimValues.landPGD = uni_values(:, 1:N_sce);
    GFLandSimValues.landProb = uni_values(:, N_sce+1:end);
    GFLandSimValues.bid_sid = [sel_BuildingInfo.building_data, build_sids];
    clear uni_values

    BuildingInfo_quick.GSSimValues = GSSimValues;
    BuildingInfo_quick.GFLatSimValues = GFLatSimValues;
    BuildingInfo_quick.GFSetSimValues = GFSetSimValues;
    BuildingInfo_quick.GFLandSimValues = GFLandSimValues;
    BuildingInfo_quick.building_edge_data = sel_BuildingInfo.building_edge_data;
    BuildingInfo_quick.wd_mean_std = sel_BuildingInfo.wd_mean_std;
    clear GSSimValues GSSimValues GFSetSimValues GFLandSimValues sel_BuildingInfo
    % bridge 
    is_bri = arrayfun(@(s) strcmp(s.Bridge, 'T'), RoadSystem.Edge);
    if size(is_bri,1)==1
        is_bri = is_bri';
    end
    bridge_edge_id = [RoadSystem.Edge(is_bri).ID]';
    zoneCells = { RoadSystem.Edge(is_bri).Zone };
    bri_zone_gmf = zeros(length(zoneCells), N_sce);
    for i = 1:length(zoneCells)
        if ~isempty(zoneCells{i})
            is_zids = ismember(SA1p0(:,1),zoneCells{i}(:,1));
            bri_zone_gmf(i,:) = max(SA1p0(is_zids, 4:end), [], 1);
        end
    end
    [uni_values, ~, bri_sids] = unique(single(bri_zone_gmf), 'rows');
    GSRoadSimValues.SA1p0 = uni_values(:, 1:N_sce);
    GSRoadSimValues.bri_sids = [bridge_edge_id, bri_sids];
    GSRoadSimValues.is_bri = is_bri;
    clear bri_zone_gmf uni_values

    % road in GF

    is_e_maj = [RoadSystem.Edge.Highway]'<= 4;
    is_bridge = arrayfun(@(s) strcmp(s.Bridge, 'T'), RoadSystem.Edge);
    if size(is_bridge,1)==1
        is_bridge = is_bridge';
    end
    is_tunnel = arrayfun(@(s) strcmp(s.Tunnel, 'T'), RoadSystem.Edge);
    if size(is_tunnel,1)==1
        is_tunnel = is_tunnel';
    end
    edge_id_type = [RoadSystem.Edge.ID]';
    edge_id_type(:,2) = 1*(is_e_maj&(~is_bridge)&(~is_tunnel)) + 2*((~is_e_maj)&(~is_bridge)&(~is_tunnel)) + 3*is_bridge+ 4*is_tunnel;
    % edge_id_type(:,3) = [RoadSystem.Edge.aoi_id]';
    % 1: major road; 2: urban road; 3: bridge; 4: tunnel
    buff = latPGD(:,4:end).*liquProb(:,4:end);
    max_liq_pgd = max(cat(3, buff, setPGD(:,4:end)), [], 3); % Maximum PGD for each road-zone and scenario
    clear latPGD setPGD 

    zoneCells = { RoadSystem.Edge.Zone };

    nZ   = numel(zoneCells);
    ids  = SA1p0(:,1);
    clear SA1p0 RoadSystem
    N_ids  = numel(ids);
    
    % 1) 估计非零数 ≈ 2 * zoneCells 长度
    nnzEst = 2 * nZ;
    
    % 2) 预分配 rowIdx 和 colIdx
    rowIdx = zeros(nnzEst,1);
    colIdx = zeros(nnzEst,1);
    
    % 3) 用一个指针 k 来填充
    k = 1;
    for i = 1:nZ
        if ~isempty(zoneCells{i})
            [~, loc] = ismember(zoneCells{i}(:,1), ids);
            loc = loc(loc>0);
            nloc = numel(loc);
            if nloc==0, continue; end
            
            % 防止超出预分配：如超出，可双倍扩容
            if k + nloc - 1 > nnzEst
                % 扩容一倍
                extend = nnzEst;
                rowIdx = [rowIdx; zeros(extend,1)];
                colIdx = [colIdx; zeros(extend,1)];
                nnzEst = nnzEst + extend;
            end
            
            rowIdx(k : k+nloc-1) = loc;
            colIdx(k : k+nloc-1) = i;
            k = k + nloc;
        end
    end
    
    % 4) 截断多余空间
    rowIdx = rowIdx(1:k-1);
    colIdx = colIdx(1:k-1);
    
    % 5) 一次性构造稀疏矩阵
    maskC = sparse(rowIdx, colIdx, true, N_ids, nZ);
    
    % （如需要逻辑型）
    maskC = logical(maskC);
    
    EdgeLiqPGD  = zeros(nZ, N_sce, 'single');
    EdgeLiqProb = zeros(nZ, N_sce, 'single');
    
    for i = 1:nZ
        m = maskC(:,i);
        if any(m)
            EdgeLiqPGD(i,:)  = max(max_liq_pgd(m,   :), [], 1);
            EdgeLiqProb(i,:)= max(liquProb(m, 4:end), [], 1);
        end
    end
    clear max_liq_pgd liquProb

    try
        [uni_values, ~, eid_sids] = unique(single([EdgeLiqPGD,EdgeLiqProb,edge_id_type(:,2)]), 'rows');
        clear EdgeLiqPGD EdgeLiqProb
        LiqRoadSimValues.liquPGD = uni_values(:, 1:N_sce);
        LiqRoadSimValues.liquProb = uni_values(:, N_sce+1:end-1);
        LiqRoadSimValues.edge_type = uni_values(:, end);
        LiqRoadSimValues.edge_id_type_sids = [edge_id_type, eid_sids];
        clear EdgeLiqPGD EdgeLiqProb uni_values
    catch ME
        if isa(EdgeLiqPGD, 'gpuArray')
            EdgeLiqPGD = gather(EdgeLiqPGD);
            EdgeLiqProb = gather(EdgeLiqProb);
        end
        [uni_values, ~, eid_sids] = unique(single([EdgeLiqPGD,EdgeLiqProb,edge_id_type(:,2)]), 'rows');
        clear EdgeLiqPGD EdgeLiqProb
        LiqRoadSimValues.liquPGD = uni_values(:, 1:N_sce);
        LiqRoadSimValues.liquProb = uni_values(:, N_sce+1:end-1);
        LiqRoadSimValues.edge_type = uni_values(:, end);
        LiqRoadSimValues.edge_id_type_sids = [edge_id_type, eid_sids];
        clear EdgeLiqPGD EdgeLiqProb uni_values
    end

    EdgeLandPGD  = zeros(nZ, N_sce, 'single');
    EdgeLandProb = zeros(nZ, N_sce, 'single');
    for i = 1:nZ
        m = maskC(:,i);
        if ~any(m), continue; end
        EdgeLandPGD(i,:)  = max(landPGD(m,   4:end), [], 1);
        EdgeLandProb(i,:) = max(landProb(m, 4:end), [], 1);
    end
    clear landPGD landProb maskC

    try
        [uni_values, ~, eid_sids] = unique(single([EdgeLandPGD,EdgeLandProb,edge_id_type(:,2)]), 'rows');
        clear EdgeLandPGD EdgeLandProb
        LandRoadSimValues.landPGD = uni_values(:, 1:N_sce);
        LandRoadSimValues.landProb = uni_values(:, N_sce+1:end-1);
        LandRoadSimValues.edge_type = uni_values(:, end);
        LandRoadSimValues.edge_id_type_sids = [edge_id_type, eid_sids];
        clear EdgeLandPGD EdgeLandProb uni_values
    catch ME
        if isa(EdgeLandPGD, 'gpuArray')
            EdgeLandPGD = gather(EdgeLandPGD);
            EdgeLandProb = gather(EdgeLandProb);
        end
        [uni_values, ~, eid_sids] = unique(single([EdgeLandPGD,EdgeLandProb,edge_id_type(:,2)]), 'rows');
        clear EdgeLandPGD EdgeLandProb
        LandRoadSimValues.landPGD = uni_values(:, 1:N_sce);
        LandRoadSimValues.landProb = uni_values(:, N_sce+1:end-1);
        LandRoadSimValues.edge_type = uni_values(:, end);
        LandRoadSimValues.edge_id_type_sids = [edge_id_type, eid_sids];
        clear EdgeLandPGD EdgeLandProb uni_values
    end

    RoadSystem_quick.GSRoadSimValues = GSRoadSimValues;
    RoadSystem_quick.LiqRoadSimValues = LiqRoadSimValues;
    RoadSystem_quick.LandRoadSimValues = LandRoadSimValues;
    RoadSystem_quick.edge_id_type = edge_id_type;
    clear GSRoadSimValues LiqRoadSimValues LandRoadSimValues
end