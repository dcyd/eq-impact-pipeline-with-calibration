function [gid1_unacc,build_unacc,gid1_unacc_init,build_unacc_init] = compute_unacc_injuries_pop(building_state, all_injure_rate, post_ptime, init_ptime,build_node_id,building_pop, batch_size)


    if nargin < 7
        batch_size = 50;  % default batch size
    end

    if ~isempty(post_ptime)
        is_post = true;
    else
        is_post = false;
    end

    if ~isempty(init_ptime)
        is_init = true;
    else
        is_init = false;
    end


    N_sce = size(building_state,2);  % number of scenarios
    N_B = size(building_state,1);  % number of buildings
    nbatch = ceil(N_sce / batch_size);

    gid1_unacc.injury_10 = zeros(1,N_sce,'single');
    gid1_unacc.injury_30 = zeros(1,N_sce,'single');
    gid1_unacc.injury_60 = zeros(1,N_sce,'single');
    gid1_unacc.injury_all = zeros(1,N_sce,'single');

    gid1_unacc.pop_10 = zeros(1,N_sce,'single');
    gid1_unacc.pop_30 = zeros(1,N_sce,'single');
    gid1_unacc.pop_60 = zeros(1,N_sce,'single');

    build_unacc.injury_10_sum = zeros(N_B,1,'single');
    build_unacc.injury_30_sum = zeros(N_B,1,'single');
    build_unacc.injury_60_sum = zeros(N_B,1,'single');
    build_unacc.pop_10 = zeros(N_B,1,'single');
    build_unacc.pop_30 = zeros(N_B,1,'single');
    build_unacc.pop_60 = zeros(N_B,1,'single');

    gid1_unacc_init.injury_10 = zeros(1,N_sce,'single');
    gid1_unacc_init.injury_30 = zeros(1,N_sce,'single');
    gid1_unacc_init.injury_60 = zeros(1,N_sce,'single');
    gid1_unacc_init.pop_10 = zeros(1,1,'single');
    gid1_unacc_init.pop_30 = zeros(1,1,'single');
    gid1_unacc_init.pop_60 = zeros(1,1,'single');

    build_unacc_init.injury_10_sum = zeros(N_B,1,'single');
    build_unacc_init.injury_30_sum = zeros(N_B,1,'single');
    build_unacc_init.injury_60_sum = zeros(N_B,1,'single');
    build_unacc_init.pop_10 = zeros(N_B,1,'single');
    build_unacc_init.pop_30 = zeros(N_B,1,'single');
    build_unacc_init.pop_60 = zeros(N_B,1,'single');

    if is_post
        build_unacc.injury_sum = zeros(N_B,1,'single');
    elseif is_init
        build_unacc_init.injury_sum = zeros(N_B,1,'single');
    end
    for b = 1:nbatch
        idx = (b-1)*batch_size + 1 : min(b*batch_size,N_sce);
        % 0: not damaged, 1: slight, 2: moderate, 3: extensive, 4: complete no collapse, 5: collapsed
        building_injuries =...
            ((building_state(:,idx,:)==1).*all_injure_rate(:,1)+... % slight
            (building_state(:,idx,:)==2).*all_injure_rate(:,2)+...  % moderate
            (building_state(:,idx,:)==3).*all_injure_rate(:,3)+... % extensive
            (building_state(:,idx,:)==4).*all_injure_rate(:,4)+... % complete no collapse
            (building_state(:,idx,:)==5).*all_injure_rate(:,5)...  % collapsed
            ).*building_pop;

        if is_post
            is_exceed = post_ptime(build_node_id(:,end),idx,:)>10;
            idx_unacc_injuries = mean(is_exceed.*building_injuries,3);
            gid1_unacc.injury_10(1,idx) = sum(idx_unacc_injuries);
            build_unacc.injury_10_sum = build_unacc.injury_10_sum + sum(idx_unacc_injuries,2);
            idx_unacc_pop = mean(is_exceed,3).*building_pop;
            gid1_unacc.pop_10(1,idx) = sum(idx_unacc_pop);
            build_unacc.pop_10 = build_unacc.pop_10 + sum(idx_unacc_pop,2);
    
            is_exceed = post_ptime(build_node_id(:,end),idx,:)>30;
            idx_unacc_injuries = mean(is_exceed.*building_injuries,3);
            gid1_unacc.injury_30(1,idx) = sum(idx_unacc_injuries);
            build_unacc.injury_30_sum = build_unacc.injury_30_sum + sum(idx_unacc_injuries,2);
            idx_unacc_pop = mean(is_exceed,3).*building_pop;
            gid1_unacc.pop_30(1,idx) = sum(idx_unacc_pop);
            build_unacc.pop_30 = build_unacc.pop_30 + sum(idx_unacc_pop,2);
    
            is_exceed = post_ptime(build_node_id(:,end),idx,:)>60;
            idx_unacc_injuries = mean(is_exceed.*building_injuries,3);
            gid1_unacc.injury_60(1,idx) = sum(idx_unacc_injuries);
            build_unacc.injury_60_sum = build_unacc.injury_60_sum + sum(idx_unacc_injuries,2);
            idx_unacc_pop = mean(is_exceed,3).*building_pop;
            gid1_unacc.pop_60(1,idx) = sum(idx_unacc_pop);
            build_unacc.pop_60 = build_unacc.pop_60 + sum(idx_unacc_pop,2);

            gid1_unacc.injury_all(1,idx) = sum(mean(building_injuries,3));
            build_unacc.injury_sum = build_unacc.injury_sum + sum(mean(building_injuries,3),2);
        end

        if is_init
            % initial accessibility
            is_exceed = init_ptime(build_node_id(:,end),1)>10;
            idx_unacc_injuries = mean(is_exceed.*building_injuries,3);
            gid1_unacc_init.injury_10(1,idx) = sum(idx_unacc_injuries);
            build_unacc_init.injury_10_sum = build_unacc_init.injury_10_sum + sum(idx_unacc_injuries,2);

            is_exceed = init_ptime(build_node_id(:,end),1)>30;
            idx_unacc_injuries = mean(is_exceed.*building_injuries,3);
            gid1_unacc_init.injury_30(1,idx) = sum(idx_unacc_injuries);
            build_unacc_init.injury_30_sum = build_unacc_init.injury_30_sum + sum(idx_unacc_injuries,2);

            is_exceed = init_ptime(build_node_id(:,end),1)>60;
            idx_unacc_injuries = mean(is_exceed.*building_injuries,3);
            gid1_unacc_init.injury_60(1,idx) = sum(idx_unacc_injuries);
            build_unacc_init.injury_60_sum = build_unacc_init.injury_60_sum + sum(building_injuries,2);

            if ~is_post
                build_unacc_init.injury_sum = build_unacc_init.injury_sum + sum(building_injuries,2);
            end

        end

        clear building_injuries is_exceed idx_unacc_injuries idx_unacc_pop
    end

    if is_init
        % initial accessibility
        is_exceed = init_ptime(build_node_id(:,end),1)>10;
        idx_unacc_pop = is_exceed.*building_pop;
        gid1_unacc_init.pop_10 = sum(idx_unacc_pop);
        build_unacc_init.pop_10 = idx_unacc_pop;

        is_exceed = init_ptime(build_node_id(:,end),1)>30;
        idx_unacc_pop = is_exceed.*building_pop;
        gid1_unacc_init.pop_30 = sum(idx_unacc_pop);
        build_unacc_init.pop_30 = idx_unacc_pop;

        is_exceed = init_ptime(build_node_id(:,end),1)>60;
        idx_unacc_pop = is_exceed.*building_pop;
        gid1_unacc_init.pop_60 = sum(idx_unacc_pop);
        build_unacc_init.pop_60 = idx_unacc_pop;
    end

end
