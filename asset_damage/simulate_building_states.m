function building_state = simulate_building_states(prob_path,rate_cd2collap,batch_size,N_sim)

    if nargin < 3
        batch_size = 50;  % default batch size
        N_sim = 10;  % default number of simulations
    end
    if nargin < 4
        N_sim = 10;  % default number of simulations
    end
    dam_states = {'slig','mode','exte','comp'};
    % load the building probability
    load(strcat(prob_path,'comp.mat'),'buildingGSfp');
    uni_building_GSfp = buildingGSfp.unique_GSfp;
    gs_bid = buildingGSfp.gs_bid;
    clear buildingGSfp

    % ---- size ----
    N_B   = size(gs_bid,1);  % number of buildings
    N_sce = size(uni_building_GSfp,2);  % number of scenarios
    nbatch = ceil(N_sce / batch_size);

    load(strcat(prob_path,'comp.mat'),'buildingGFfp');
    isGFstruct = isstruct(buildingGFfp);
    if ~isGFstruct
        buildingGFfp = single(buildingGFfp);            % 若已是矩阵
    end
    % initialize
    building_state = zeros(N_B,N_sce,N_sim,'uint8');  % damage states
    % 0: not damaged, 1: slight, 2: moderate, 3: extensive, 4: complete no collapse, 5: collapsed
    % generate random numbers
    pro_build_rand = rand(N_B, batch_size, N_sim,'single');
    for b = 1:nbatch
        idx = (b-1)*batch_size + 1 : min(b*batch_size,N_sce);
        i_nsce = length(idx);

        GSb = uni_building_GSfp(gs_bid,idx);
        if isGFstruct
            GFb = single(buildingGFfp.zone_GFfp_cp(buildingGFfp.bid_zid, idx));
        else
            GFb = buildingGFfp(:,idx);
        end

        mix_comp_prob = GSb + GFb - GSb.*GFb;
        mix_comp_prob_collapsed = gather(mix_comp_prob.*rate_cd2collap);  % complete to collapsed ratio  
        
        for si = 1:N_sim
            building_state(:,idx,si) = uint8(pro_build_rand(:,1:i_nsce,si)<mix_comp_prob_collapsed(:,:));
            building_state(:,idx,si) = building_state(:,idx,si)+uint8(pro_build_rand(:,1:i_nsce,si)<mix_comp_prob(:,:));
        end
    end
    clear uni_building_GSfp buildingGFfp GSb GFb mix_comp_prob mix_comp_prob_collapsed rate_cd2collap gs_bid

    for dsi = 3:-1:1
        load(strcat(prob_path,dam_states{dsi},'.mat'), ...
                'buildingGSfp');
        uni_building_GSfp = buildingGSfp.unique_GSfp;
        gs_bid = buildingGSfp.gs_bid;
        clear buildingGSfp

        load(strcat(prob_path,dam_states{dsi},'.mat'), ...
                'buildingGFfp');
        if ~isGFstruct, buildingGFfp = single(buildingGFfp); end

        for b = 1:nbatch
            idx = (b-1)*batch_size + 1 : min(b*batch_size,N_sce);
            i_nsce = length(idx);
            GS_new = uni_building_GSfp(gs_bid,idx);

            if ~isempty(buildingGFfp)
                if isGFstruct
                    if isfield(buildingGFfp,"zone_GFfp_cp")
                        GF_new = single(buildingGFfp.zone_GFfp_cp( ...
                                        buildingGFfp.bid_zid, idx));
                    elseif isfield(buildingGFfp,"zone_GFfp_ex")
                        GF_new = single(buildingGFfp.zone_GFfp_ex( ...
                                        buildingGFfp.bid_zid, idx));
                    end
                else
                    GF_new = buildingGFfp(:,idx);    
                end
                % ---- GS+GF 增量 ----
                mix_new = GS_new + GF_new - GS_new.*GF_new;
                clear GF_new
            else
                mix_new = GS_new;
            end
            for si = 1:N_sim
                building_state(:,idx,si) = building_state(:,idx,si)+uint8(pro_build_rand(:,1:i_nsce,si)<mix_new(:,:));
            end
        end
        clear uni_building_GSfp GS_new mix_new
    end
end
