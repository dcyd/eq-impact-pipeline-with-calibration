function demo_accessibility()
    % DEMO_LOSS_FATALITY  Quick-start example for the seismic impact pipeline.
    % This demo is set in Naypyitaw, Myanmar (the GID1 region in GADM data), in the 2025 earthquake.

    fprintf('Quick-start demo: accessibility for Naypyitaw, Myanmar in the 2025 earthquake\n\n');
    %% 0. Environment & paths

    % Check for GPU availability
    try
        if gpuDeviceCount > 0
            try
                gpuDevice(1); % Try to select the first GPU to confirm availability
                isGpu = true;
            catch
                isGpu = false;
            end
        else
            isGpu = false;
        end
    catch
        isGpu = false;
    end

    % file path
    demo_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(demo_dir, 'asset_damage'));
    addpath(fullfile(demo_dir, 'socio_economic_impact'));
    addpath(fullfile(demo_dir, 'utils'));
    addpath(fullfile(demo_dir, 'utils\matlab_bgl'))
    data_dir = fullfile(demo_dir, 'data');
    buildingInfoFile       = fullfile(data_dir, 'BuildingInfo_ID10.mat');
    seismicScenarioFile    = fullfile(data_dir, 'SeismicScenario_ID10.mat');
    buildQuickSimData_file = fullfile(data_dir, 'BuildQuickSimData_ID10.mat');
    buildProDataDefault_file      = fullfile(data_dir, 'BuildProb_ID10_default_');
    buildProDataCorrect_file      = fullfile(data_dir, 'BuildProb_ID10_correct_');

    roadSystemFile = fullfile(data_dir,"RoadSystem_ID10.mat");
    roadQuickSimData_file = fullfile(data_dir, 'RoadQuickSimData_ID10.mat');
    roadProDataDefault_file      = fullfile(data_dir, 'RoadProb_ID10_default.mat');
    roadProDataCorrect_file      = fullfile(data_dir, 'RoadProb_ID10_correct.mat');
    roadMatrix_file = fullfile(data_dir, 'RoadSystem_matrix_ID10.mat');

    selBuildStateDefault_file = fullfile(data_dir, "selstate_ID10_default.mat");
    TTDefault_file = fullfile(data_dir,"TT_ID10_default.mat");
    BuildStateDefault_file = fullfile(data_dir, "BuildingState_ID10_default.mat");

    selBuildStateCorrect_file = fullfile(data_dir, "selstate_ID10_correct.mat");
    TTCorrect_file = fullfile(data_dir,"TT_ID10_correct.mat");
    BuildStateCorrect_file = fullfile(data_dir, "BuildingState_ID10_correct.mat");

    hazusRates_file = fullfile(data_dir,'hazus_rates_pkg.mat');
    buildingHazusIDX_file = fullfile(data_dir,'building_hazus_casualty_idx_ID10.mat');

    % correction factors
    correctFactors = [1.43, 1.01, 1.03, 1.66 1.52 2.46 1.27 1.30 1.00 1.08 1.01];  % calibrated
    defaultFactors = ones(1,11);  % default

    %% 1. Load inputs

    if ~exist(buildQuickSimData_file, "file")
        if ~exist('BuildingInfo', 'var')
            load(buildingInfoFile, 'BuildingInfo');
        end
        if ~exist('SeismicScenario', 'var')
            load(seismicScenarioFile, 'SeismicScenario');
        end

        [BuildingInfo_quick] = generate_building_quick_sim_data(BuildingInfo,SeismicScenario,[],isGpu);
        save(buildQuickSimData_file, "BuildingInfo_quick","-v7.3");
    end
    
    if ~ exist(roadQuickSimData_file,"file")    
        if ~exist('BuildingInfo', 'var')
            load(buildingInfoFile, 'BuildingInfo');
        end
        if ~exist('SeismicScenario', 'var')
            load(seismicScenarioFile, 'SeismicScenario');
        end
        if ~exist('RoadSystem', 'var')
            load(roadSystemFile, "RoadSystem");
        end

        [selBuildingInfo_quick,RoadSystem_quick] = generate_road_quick_sim_data(BuildingInfo, RoadSystem, SeismicScenario,isGpu);
        save(roadQuickSimData_file, "selBuildingInfo_quick", "RoadSystem_quick","-v7.3");
    end
    clear BuildingInfo SeismicScenario RoadSystem

    load(buildQuickSimData_file, 'BuildingInfo_quick');
    load(roadQuickSimData_file, "selBuildingInfo_quick", "RoadSystem_quick");
    load(fullfile(data_dir, 'building_pop_lethality_ID10.mat'),'building_pop_lethality');
    load(fullfile(data_dir,'Build_GS_Fra_ID10.mat'),'fra_params');

    %% 2. Compute building-level damage probabilities in default pipeline

    if ~exist(strcat(buildProDataDefault_file, 'comp.mat'),'file')
        FreeParams.cfBuildingLatMedia = defaultFactors(1);
        FreeParams.cfBuildingSetMedia = defaultFactors(2);
        FreeParams.cfBuildingLandMedia = defaultFactors(3);
        FreeParams.cfBuildGSMedia = defaultFactors(4);
    
        dam_states = {"slig","mode","exte","comp"};
        for dsi = 1:4
            ds = dam_states{dsi};
            FreeParams.state = ds;
            ds_BuildingInfo_quick = BuildingInfo_quick;
            ds_BuildingInfo_quick.GSSimValues.fra_media = ds_BuildingInfo_quick.GSSimValues.fra_media(:,dsi);
    
            [buildingGSfp,buildingGFfp] = calculate_building_damage_prob(ds_BuildingInfo_quick, FreeParams);
            save(strcat(buildProDataDefault_file, ds, '.mat'), "buildingGSfp","buildingGFfp","-v7.3");
    
            clear buildingGSfp buildingGFfp ds_BuildingInfo_quick
        end
    end

    %% 3. Compute building-level damage probabilities given correction factors
    if ~exist(strcat(buildProDataCorrect_file, 'comp.mat'),'file')
        FreeParams.cfBuildingLatMedia = correctFactors(1);
        FreeParams.cfBuildingSetMedia = correctFactors(2);
        FreeParams.cfBuildingLandMedia = correctFactors(3);
        FreeParams.cfBuildGSMedia = correctFactors(4);
    
        dam_states = {"slig","mode","exte","comp"};
        for dsi = 1:4
            ds = dam_states{dsi};
            FreeParams.state = ds;
            ds_BuildingInfo_quick = BuildingInfo_quick;
            ds_BuildingInfo_quick.GSSimValues.fra_media = ds_BuildingInfo_quick.GSSimValues.fra_media(:,dsi);
            try
                [buildingGSfp,buildingGFfp] = calculate_building_damage_prob(ds_BuildingInfo_quick, FreeParams);
                save(strcat(buildProDataCorrect_file, ds, '.mat'), "buildingGSfp","buildingGFfp","-v7.3");
            catch
                sprintf("error of the city %d in the %s state for the baseline results",GID_code,ds)
            end
            clear buildingGSfp buildingGFfp ds_BuildingInfo_quick
        end
    end

    %% 4. Generate road-segment-level damage probability data - default

    if ~exist(roadProDataDefault_file,"file")
        FreeParams.cfBuildingLatMedia = defaultFactors(1);
        FreeParams.cfBuildingSetMedia = defaultFactors(2);
        FreeParams.cfBuildingLandMedia = defaultFactors(3);
        FreeParams.cfBuildGSMedia = defaultFactors(4);
        FreeParams.cfBridgeGSMedia = defaultFactors(5);
        FreeParams.cfUrbanRoadGFMedia = defaultFactors(6);
        FreeParams.cfMajorRoadGFMedia = defaultFactors(7);
        FreeParams.cfBridgeGFMedia = defaultFactors(8);
        FreeParams.cfTunnelGFMedia = defaultFactors(9);
        FreeParams.cfWDMean = defaultFactors(10);
        FreeParams.cfWDSTD = defaultFactors(11);

        [bridge_GSfp,edge_Block_GSfp,edge_GFfp,edge_Block_GFfp] = calculate_road_damage_prob(selBuildingInfo_quick, RoadSystem_quick, FreeParams);

        save(roadProDataDefault_file, "bridge_GSfp","edge_Block_GSfp","edge_GFfp","edge_Block_GFfp","-v7.3");
        clear bridge_GSfp edge_Block_GSfp edge_GFfp edge_Block_GFfp
    end

    %% 5. Generate road-segment-level damage probability data - calibrated

    if ~exist(roadProDataCorrect_file,"file")
        FreeParams.cfBuildingLatMedia = correctFactors(1);
        FreeParams.cfBuildingSetMedia = correctFactors(2);
        FreeParams.cfBuildingLandMedia = correctFactors(3);
        FreeParams.cfBuildGSMedia = correctFactors(4);
        FreeParams.cfBridgeGSMedia = correctFactors(5);
        FreeParams.cfUrbanRoadGFMedia = correctFactors(6);
        FreeParams.cfMajorRoadGFMedia = correctFactors(7);
        FreeParams.cfBridgeGFMedia = correctFactors(8);
        FreeParams.cfTunnelGFMedia = correctFactors(9);
        FreeParams.cfWDMean = correctFactors(10);
        FreeParams.cfWDSTD = correctFactors(11);

        [bridge_GSfp,edge_Block_GSfp,edge_GFfp,edge_Block_GFfp] = calculate_road_damage_prob(selBuildingInfo_quick, RoadSystem_quick, FreeParams);

        save(roadProDataCorrect_file, "bridge_GSfp","edge_Block_GSfp","edge_GFfp","edge_Block_GFfp","-v7.3");
        clear bridge_GSfp edge_Block_GSfp edge_GFfp edge_Block_GFfp
    end


    %% 6. Simulate building state and travel time to hospitals - default
    rate_cd2collap = fra_params.rate_cd2collap; clear fra_params

    load(roadMatrix_file,"RoadSystem_mat")

    prob_path = buildProDataDefault_file;


    if ~exist(BuildStateDefault_file,"file")
        building_state = simulate_building_states(prob_path,rate_cd2collap);
        save(BuildStateDefault_file, "building_state", "-v7.3");
    else
        load(BuildStateDefault_file, "building_state");
    end


    if ~exist(TTDefault_file,"file")

        % building damage state simulation
        if ~exist(selBuildStateDefault_file,"file")
            [~,sel_build_tag] = ismember(selBuildingInfo_quick.building_edge_data(:,1),building_pop_lethality(:,1));
            sel_build_state.sel_edge_build_state = building_state(sel_build_tag,:,:);

            hos_node_build = identify_hos_node_within_bound(RoadSystem_mat.hos_node(:,2:3),building_pop_lethality(:,1:3));
            [~,hos_build_ids] = ismember(hos_node_build(:,end),building_pop_lethality(:,1));  % building idx
            sel_build_state.sel_hos_build_state = building_state(hos_build_ids,:,:);
            clear building_state
            save(selBuildStateDefault_file, "sel_build_state", "-v7.3");
        else
            load(selBuildStateDefault_file, "sel_build_state");
        end

        % road blockage probability and simulations
        % ---- initialize the road-building probability of building debris blockage----

        edge_ids = RoadSystem_quick.edge_id_type(:,1); clear RoadSystem_quick
    
        load(roadProDataDefault_file, "bridge_GSfp","edge_GFfp");
        road_prob.bridge_GSfp = bridge_GSfp;
        road_prob.edge_GFfp = edge_GFfp;
        clear bridge_GSfp edge_GFfp
    
        FreeParams.cfWDMean = defaultFactors(10);
        FreeParams.cfWDSTD = defaultFactors(11);

        [post_ptime,init_ptime] = ...
            simulate_tt_given_building_state(selBuildingInfo_quick,edge_ids,...
            RoadSystem_mat,sel_build_state,road_prob,FreeParams);
        
        save(TTDefault_file, "post_ptime", "init_ptime", "-v7.3");
        clear init_ptime post_ptime
    end

    %% 7. Simulate building state and travel time to hospitals - calibrated
    prob_path = buildProDataCorrect_file;

    if ~exist(BuildStateCorrect_file,"file")
        building_state = simulate_building_states(prob_path,rate_cd2collap);
        save(BuildStateCorrect_file, "building_state", "-v7.3");
    else
        load(BuildStateCorrect_file, "building_state");
    end

    if ~exist(TTCorrect_file,"file")

        % building damage state simulation
        if ~exist(selBuildStateCorrect_file,"file")
            [~,sel_build_tag] = ismember(selBuildingInfo_quick.building_edge_data(:,1),building_pop_lethality(:,1));
            sel_build_state.sel_edge_build_state = building_state(sel_build_tag,:,:);

            hos_node_build = identify_hos_node_within_bound(RoadSystem_mat.hos_node(:,2:3),building_pop_lethality(:,1:3));
            [~,hos_build_ids] = ismember(hos_node_build(:,end),building_pop_lethality(:,1));  % building idx
            sel_build_state.sel_hos_build_state = building_state(hos_build_ids,:,:);
            clear building_state
            save(selBuildStateCorrect_file, "sel_build_state", "-v7.3");
        else
            load(selBuildStateCorrect_file, "sel_build_state");
        end

        % road blockage probability and simulations
        % ---- initialize the road-building probability of building debris blockage----

        load(roadProDataCorrect_file, "bridge_GSfp","edge_GFfp");
        road_prob.bridge_GSfp = bridge_GSfp;
        road_prob.edge_GFfp = edge_GFfp;
        clear bridge_GSfp edge_GFfp
    
        FreeParams.cfWDMean = correctFactors(10);
        FreeParams.cfWDSTD = correctFactors(11);

        [post_ptime,init_ptime] = ...
            simulate_tt_given_building_state(selBuildingInfo_quick,edge_ids,...
            RoadSystem_mat,sel_build_state,road_prob,FreeParams);
        
        save(TTCorrect_file, "post_ptime", "init_ptime", "-v7.3");
        clear init_ptime post_ptime
    end

    %% 8. Simulate accessibility results - default

    L = load(hazusRates_file,'R','rateCols');
    I = load(buildingHazusIDX_file,'hazus_idx');
    injure_rates = expand_rates(I.hazus_idx, L.R);
    injure_rates = injure_rates/100;
    all_injure_rate = single([ ...
            sum(injure_rates(:,1:3),  2), ... % 1: slight
            sum(injure_rates(:,5:7),  2), ... % 2: moderate
            sum(injure_rates(:,9:11), 2), ... % 3: extensive
            sum(injure_rates(:,13:15),2), ... % 4: complete no collapse
            sum(injure_rates(:,17:19),2)  ... % 5: collapsed
        ]);
    clear injure_rates I L
    build_node_id = RoadSystem_mat.build_node_id;

    load(BuildStateDefault_file, "building_state");
    load(TTDefault_file, "post_ptime","init_ptime");

    [gid1Unacc_default, ~, ~, ~] = compute_unacc_injuries_pop(building_state, all_injure_rate, ...
            post_ptime, init_ptime, build_node_id, building_pop_lethality(:,7));
    
    clear post_ptime init_ptime building_state

    %% 9. Simulate accessibility results - calibrated

    load(BuildStateCorrect_file, "building_state");
    load(TTCorrect_file, "post_ptime","init_ptime");

    [gid1Unacc_correct, ~, ~, ~] = compute_unacc_injuries_pop(building_state, all_injure_rate, ...
            post_ptime, init_ptime, build_node_id, building_pop_lethality(:,7));
    
    clear post_ptime init_ptime building_state all_injure_rate
  
    %% 10. Portfolio summary
    fprintf('\n=== Portfolio statistics ===\n');

    fprintf('Total injuries with travel time to hospital more than 60 minutes in default pipeline (median): %.2e, 95%% range: [%.2e, %.2e]\n', ...
        median(gid1Unacc_default.injury_60), prctile(gid1Unacc_default.injury_60, 2.5), prctile(gid1Unacc_default.injury_60, 97.5));

    fprintf('Total injuries in default pipeline (median): %.2e, 95%% range: [%.2e, %.2e]\n', ...
        median(gid1Unacc_default.injury_all), prctile(gid1Unacc_default.injury_all, 2.5), prctile(gid1Unacc_default.injury_all, 97.5));

    fprintf('Total injuries with travel time to hospital more than 60 minutes in default pipeline (median): %.2e, 95%% range: [%.2e, %.2e]\n', ...
        median(gid1Unacc_correct.injury_60), prctile(gid1Unacc_correct.injury_60, 2.5), prctile(gid1Unacc_correct.injury_60, 97.5));

    fprintf('Total injuries in default pipeline (median): %.2e, 95%% range: [%.2e, %.2e]\n', ...
        median(gid1Unacc_correct.injury_all), prctile(gid1Unacc_correct.injury_all, 2.5), prctile(gid1Unacc_correct.injury_all, 97.5));

end