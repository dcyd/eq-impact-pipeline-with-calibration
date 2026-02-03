function demo_loss_fatality()
    % DEMO_LOSS_FATALITY  Quick-start example for the seismic impact pipeline.
    % This demo is set in Naypyitaw, Myanmar (the GID1 region in GADM data), in the 2025 earthquake.

    % 0. Environment & paths
    fprintf('Quick-start demo: direct economic loss and fatalities for Naypyitaw, Myanmar in the 2025 earthquake\n\n');
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

    demo_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(demo_dir, 'asset_damage'));
    addpath(fullfile(demo_dir, 'socio_economic_impact'));

    data_dir = fullfile(demo_dir, 'data');
    buildingInfoFile       = fullfile(data_dir, 'extended_20_BuildingInfo_ID10.mat');
    seismicScenarioFile    = fullfile(data_dir, 'extended_20_SeismicScenario_ID10.mat');
    buildQuickSimData_file = fullfile(data_dir, 'extended_20_BuildQuickSimData_ID10.mat');
    buildProDataDefault_file      = fullfile(data_dir, '\extended_20_BuildProb_ID10_default_');
    buildProDataCorrect_file      = fullfile(data_dir, '\extended_20_BuildProb_ID10_correct_');

    load(fullfile(data_dir, 'extended_20_building_pop_lethality_ID10.mat'),'building_pop_lethality');
    load(fullfile(data_dir, 'extended_20_building_vulnerability_ID10.mat'),'building_vulnerability');

    loss_file = fullfile(data_dir, 'demo00_economic_loss.mat');
    fatality_file = fullfile(data_dir, 'demo00_fatalities.mat');

    %% 1. Load the packaged demo data
    if ~exist('BuildingInfo', 'var')
        load(buildingInfoFile, 'BuildingInfo');
    end
    if ~exist('SeismicScenario', 'var')
        load(seismicScenarioFile, 'SeismicScenario');
    end
    %% 2. Generate accelerated simulation inputs
    if ~exist(buildQuickSimData_file, "file")
        [BuildingInfo_quick] = generate_building_quick_sim_data(BuildingInfo,SeismicScenario,[],isGpu);
        save(buildQuickSimData_file, "BuildingInfo_quick","-v7.3");
    else
        load(buildQuickSimData_file, 'BuildingInfo_quick');
    end
    clear BuildingInfo SeismicScenario

    
    %% 3. Compute building-level damage probabilities in default pipeline
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


    %% 4. Compute building-level damage probabilities given correction factors
    correctFactors = [1.43, 1.01, 1.03, 1.66 1.52 2.46 1.27 1.30 1.00 1.08 1.01];  % calibrated

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


    %% 4. Propagate to direct building economic loss

    buildingStruCost = building_pop_lethality(:,22);
    buildingNonstruCost = building_pop_lethality(:,23);
    buildingContCost = building_pop_lethality(:,24);

    regionLossPerSce = struct();
    buildingLRSceSum = struct();
    % structural calibration
    [regionLossPerSce.struDefault,buildingLRSceSum.struDefault] = ...
        calculate_build_structural_loss(buildProDataDefault_file,buildingStruCost);

    [regionLossPerSce.struCorrected,buildingLRSceSum.struCorrected] = ...
        calculate_build_structural_loss(buildProDataCorrect_file,buildingStruCost);

    load(buildQuickSimData_file, 'BuildingInfo_quick');
    GSSimValues = BuildingInfo_quick.GSSimValues;
    clear BuildingInfo_quick
    [regionLossPerSce.nonStrCont,buildingLRSceSum.nonStrCont] =...
        calculate_build_other_loss(GSSimValues,building_vulnerability,[buildingNonstruCost,buildingContCost]);

    save(loss_file,'regionLossPerSce','buildingLRSceSum','-v7.3');

    %% 5. Propagate to fatalities

    load(strcat(buildProDataCorrect_file,'comp.mat'),'buildingGSfp','buildingGFfp');
    if ~exist('BuildingInfo', 'var')
        load(buildingInfoFile, 'BuildingInfo');
    end
    BuildingInfo.building_pop_lethality = building_pop_lethality;

    [regionBuildColl,regionFatality,buildFatalitySum] = ...
        calculate_collapse_fatality(buildingGSfp,buildingGFfp,BuildingInfo,50);

    save(fatality_file,'regionBuildColl','regionFatality','buildFatalitySum','-v7.3');
    clear buildingGSfp buildingGFfp BuildingInfo
    %% 6. Portfolio summary & visualization
    fprintf('\n=== Portfolio statistics ===\n');

    lossTotalDefault = regionLossPerSce.nonStrCont.lossCont + ...
        regionLossPerSce.nonStrCont.lossNonStr + ...
        regionLossPerSce.struDefault.lossGSGF; 

    post2Base = regionLossPerSce.struCorrected.lossGSGF ./...
        regionLossPerSce.struDefault.lossGSGF;

    lossTotalCorrected = regionLossPerSce.nonStrCont.lossCont.*post2Base + ...
        regionLossPerSce.nonStrCont.lossNonStr.*post2Base + ...
        regionLossPerSce.struCorrected.lossGSGF;

    fprintf('Total direct building loss (median): %.2e USD, 95%% range (loss): [%.2e, %.2e] USD\n', ...
        median(lossTotalCorrected), prctile(lossTotalCorrected, 2.5), prctile(lossTotalCorrected, 97.5));

    fprintf('Total fatalities (median): %.1f people, 95%% range (fatalities): [%.1f, %.1f] people\n', ...
        median(regionFatality), prctile(regionFatality, 2.5), prctile(regionFatality, 97.5));

    fprintf('Done. Inspect the generated figures for spatial patterns.\n');
end