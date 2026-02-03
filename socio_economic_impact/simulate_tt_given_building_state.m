function [post_ptime,init_ptime] = simulate_tt_given_building_state(BuildingInfo,edge_ids,RoadSystem,sel_build_state,road_prob,params)

    if isfield(params, 'cfWDMean')
        cf_wd_mean = params.cfWDMean;
    else
        cf_wd_mean = 1;
    end
    if isfield(params, 'cfWDSTD')
        cf_wd_std = params.cfWDSTD;
    else
        cf_wd_std = 1;
    end
    if isfield(params, 'thres_width')
        thres_width = params.thres_width;
    else
        thres_width = 3.5;
    end
    if isfield(params, 'batch_size')
        batch_size = params.batch_size;
    else
        batch_size = 50;
    end
    if isfield(params, 'hos_build_state')
        hos_build_state = params.hos_build_state;
    else
        hos_build_state = 4;
    end
    if isfield(params, 'hos_build_ratio')
        hos_build_ratio = params.hos_build_ratio;
    else
        hos_build_ratio = [1.0, 0.7, 0.5, 0, 0];
    end
    % {1.0, 0.9, 0.7, 0, 0} (positive) 
    % and {1.0, 0.5, 0.3, 0, 0} (Negative)
    % 0: not damaged, 1: slight, 2: moderate, 3: extensive, 4: complete no collapse, 5: collapsed

    sel_edge_build_state = sel_build_state.sel_edge_build_state;
    sel_hos_build_state = sel_build_state.sel_hos_build_state;
    clear sel_build_state
    N_sce = size(sel_hos_build_state, 2);
    N_sim = size(sel_hos_build_state, 3);
    nbatch = ceil(N_sce / batch_size);

    map_hos_build_func  = [hos_build_ratio(:); 0];                          % 0..4 from ratio; 5 (collapsed) → 0
    sel_hos_func = rand(size(sel_hos_build_state)) < map_hos_build_func(sel_hos_build_state+1);
    clear sel_hos_build_state

    building_edge_data = BuildingInfo.building_edge_data;
    wd_mean_std = BuildingInfo.wd_mean_std;
    wd_mean_std(:,1) = wd_mean_std(:,1)*cf_wd_mean;
    wd_mean_std(:,2) = wd_mean_std(:,2)*cf_wd_std;
    blockage_threshold = building_edge_data(:,3) + building_edge_data(:,6) - thres_width;
    edge_build_fp_wd = single(1-normcdf(blockage_threshold, wd_mean_std(:,1), wd_mean_std(:,2)));
    % [~,sel_build_tag] = ismember(BuildingInfo.building_edge_data(:,1),building_idloglat(:,1));
    [~,edge_tag] = ismember(building_edge_data(:, 2), edge_ids);
    % edge_tag_2_batch = repelem((1:batch_size)', numel(edge_tag));

    clear BuildingInfo building_edge_data wd_mean_std blockage_threshold

    % load road-related data
    in_bd_edge = RoadSystem.edge_data(:,7);
    N_edge = size(in_bd_edge,1);

    % initial accessibility
    N = size(RoadSystem.node_data,1);
    init_edge_time = RoadSystem.edge_data(:, 4) ./ (RoadSystem.edge_data(:, 6) * 1000 / 60);  % Initial travel time (in minutes)
    new_vector=zeros(N,1);
    new_vector(RoadSystem.hos_node(:,4))=0.0001;

    % travel time from node to hosiptal nodes before earthquake (minute)
    fnet=sparse([RoadSystem.edge_data(:,2);RoadSystem.edge_data(:,3)],[RoadSystem.edge_data(:,3);RoadSystem.edge_data(:,2)],...
        [init_edge_time;init_edge_time],N,N);
    [ptime,~]=shortest_paths([fnet new_vector;new_vector' 0],N+1);
    init_ptime = reshape(ptime(1:N),N,1)-0.0001;
    roadRand = rand(N_edge,batch_size,N_sim,'single');

    post_ptime = zeros(N,N_sce,N_sim,'uint8')+255;

    bridge_GSfp = road_prob.bridge_GSfp;
    edge_GFfp = road_prob.edge_GFfp;
    clear road_prob

    for b = 1:nbatch
        idx = (b-1)*batch_size + 1 : min(b*batch_size,N_sce);

        edge_tag_2_batch = repelem((1:length(idx))', numel(edge_tag));
        edge_tag_expanded_batch = repmat(edge_tag, 1, length(idx));

        % road blockage probability
        edge_block_prob = zeros(numel(edge_ids), length(idx), N_sim,"single");
        edge_build_block_prob = single(edge_build_fp_wd.*(sel_edge_build_state(:,idx,:)==0)); % 
        for si=1:N_sim
            si_edge_build_block_prob = edge_build_block_prob(:,:,si);
            edge_block_prob(:,:,si) = single(1 - accumarray([edge_tag_expanded_batch(:), edge_tag_2_batch(:)], ...
                                                    1 - si_edge_build_block_prob(:), ...
                                                    [numel(edge_ids), length(idx)], ...
                                                    @prod));
        end
        clear edge_tag_expanded_batch edge_build_block_prob

        edge_block_prob(~ismember((1:numel(edge_ids))', edge_tag), :) = 0;
        bridge_GSfp_all = (1-single(bridge_GSfp.unique_GSfp(bridge_GSfp.gs_bid,idx)));
        edge_block_prob = (1-single(edge_block_prob));
        bridge_and_block = bridge_GSfp_all.*edge_block_prob;
        clear edge_block_prob bridge_GSfp_all
        edge_gf = (1-single(edge_GFfp.unique_GFfp(edge_GFfp.gf_eid,idx)));
        roadFP = 1-bridge_and_block.*edge_gf;
        roadFP = 1-(1-roadFP).^(1/3);
        clear bridge_and_block edge_gf

        % simulated road blockage scenarios
        roadda = roadRand(:,1:length(idx),:) < roadFP;
        roadda(in_bd_edge==0,:) = 0;
        clear roadFP
        for sci = 1:length(idx)
            for ri=1:N_sim               
                if any(roadda(:,sci,ri))
                    edge_time = init_edge_time;
                    edge_time(roadda(:,sci,ri)) = inf;
                    fnet=sparse([RoadSystem.edge_data(:,2);RoadSystem.edge_data(:,3)],[RoadSystem.edge_data(:,3);RoadSystem.edge_data(:,2)],...
                        [edge_time;edge_time],N,N);
                    post_vector=new_vector;  
                    % 0: not damaged, 1: slight, 2: moderate, 3: extensive, 4: complete no collapse, 5: collapsed
                    % post_vector(RoadSystem.hos_node(sel_hos_build_state(:,idx(sci),ri)>=hos_build_state,4))=0; %TODO
                    post_vector(RoadSystem.hos_node(~sel_hos_func(:,idx(sci),ri),4))=0; %TODO
                    [ptime,~]=shortest_paths([fnet post_vector;post_vector' 0],N+1);
                    post_ptime(:,idx(sci),ri) = reshape(ptime(1:N),N,1)-0.0001;
                else
                    post_ptime(:,idx(sci),ri) = init_ptime;
                end                                
            end
        end
        clear edge_time fnet ptime
    end

end
