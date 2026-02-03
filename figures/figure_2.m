%% Fig. 2a and 2b are based on QGIS.

%% Fig. 2c

clear

result_dir = "figures\results";

fileB = fullfile(result_dir,'allChains_building.mat');   % columns 1–4
fileR = fullfile(result_dir,'allChains_road.mat');       % columns 1–7

labelsB = {'C_{B,lat}', 'C_{B,ver}', ...
           'C_{B,land}', 'C_{B,gm}'};
labelsR = {'C_{Br,gm}','C_{Ur,gf}','C_{Ma,gf}', ...
           'C_{Br,gf}','C_{Tu,gf}','C_{μ}','C_{σ}'};

paramNames = {'{C}_{B,lat}', ...
                  '{C}_{B,ver}', ...
                  '{C}_{B,land}', ...
                  '{C}_{B,gm}', ...
                  '{C}_{Br,gm}', ...
                  '{C}_{Ur,gf}', ...
                  '{C}_{Ma,gf}', ...
                  '{C}_{Br,gf}', ...
                  '{C}_{Tu,gf}', ...
                  '{C}_{\mu}', ...
                  '{C}_{\sigma}'};

re_order = [4,1,2,3,5,8,9,7,6,10,11];

% Plot options
halfEye      = false;         % true = half-eye; false = full violin
violinWidth  = 0.35;          % max half-thickness in x-units
alphaFill    = 0.28;          % violin fill transparency
lwInterval   = 2.2;           % line width for 95% CrI
msMedian     = 10;            % marker size for median
cBuild       = [0.20 0.45 0.75];
cRoad        = [0.95 0.60 0.10];

% Load and extract arrays
B = load(fileB);
R = load(fileR);
A_b = largestNumericMatrix(B);
A_r = largestNumericMatrix(R);

% Ensure params are columns
if size(A_b,2) < 4 && size(A_b,1) >= 4, A_b = A_b.'; end
if size(A_r,2) < 7 && size(A_r,1) >= 7, A_r = A_r.'; end

samplesB = A_b(:, 1:4);        % N_b x 4
samplesR = A_r(:, 1:7);        % N_r x 7

samplesB = samplesB(:,re_order(1:4));
samplesR = samplesR(:,re_order(5:end)-4);
samplesR(:,5) = samplesR(:,5)-1;
% Combine for plotting
labels = [labelsB, labelsR];
labels = paramNames;
labels = labels(re_order);

group  = [repmat("Building",1,numel(labelsB)), repmat("Road",1,numel(labelsR))];
S      = [num2cell(samplesB,1), num2cell(samplesR,1)];   % 1×(4+7) cell; each cell is a vector

% Plot
figure('Color','w','Position',[100 100 980 620]); hold on;

% Unity line (now along y, since posterior on y-axis)
yline(1,'k--','LineWidth',1);

% X positions (1 at left)
n    = numel(S);
xpos = 1:n;

S_pos = cell(size(S));
lo95  = nan(1,n);
hi95  = nan(1,n);
med   = nan(1,n);

for i = 1:n
    xi = S{i}(:);
    xi = xi(xi > 0);

    if isempty(xi)
        warning('Parameter %d has no positive samples; skipped in log plot.', i);
        S_pos{i} = NaN;
        continue;
    end

    S_pos{i} = xi;
    lo95(i)  = prctile(xi, 2.5);
    hi95(i)  = prctile(xi, 97.5);
    med(i)   = median(xi);
end

% Draw vertical violins
for i = 1:n
    this = S_pos{i};
    if all(isnan(this)), continue; end

    col  = cBuild; 
    if group(i) == "Road", col = cRoad; end
    drawViolinV(this, xpos(i), violinWidth, col, alphaFill, halfEye);
end

% Overlay CrIs (vertical) and medians
for i = 1:n
    if isnan(lo95(i)) || isnan(hi95(i)) || isnan(med(i)), continue; end

    line([xpos(i), xpos(i)], [lo95(i), hi95(i)], ...
        'Color',[0 0 0], 'LineWidth', lwInterval);
    plot(xpos(i), med(i), 'r_', 'MarkerSize', msMedian, 'LineWidth', 1.5);
end

% Axes & labels (parameters on x-axis)
set(gca, 'XTick', xpos, 'XTickLabel', labels, ...
         'XTickLabelRotation', 45, ...
         'Box','off');

ylabel('Posterior');

% ----- robust global y-limits
valid = ~isnan(lo95) & ~isnan(hi95);
ymin  = min(lo95(valid));
ymax  = max(hi95(valid));

logmin = log10(ymin);
logmax = log10(ymax);
dlog   = logmax - logmin;
logmin = logmin - 0.05*dlog;
logmax = logmax + 0.05*dlog;

ylim([10^logmin, 10^logmax]);

xlim([0.5, n+0.5]);

%% Figs. 2d and 2e

load(fullfile(result_dir,"AOI_build_road_post_results.mat"), 'building_aoi_sim', ...
    'road_aoi_sim', 'building_aoi_obs', 'road_aoi_obs', ...
    'building_aoi_sim_defau','road_aoi_sim_defau');
load(fullfile(result_dir,"aoi_info.mat"));

plot_aoi_interval_compare_defau(building_aoi_sim, building_aoi_obs, building_aoi_sim_defau,...
                          road_aoi_sim, road_aoi_obs, road_aoi_sim_defau, aoi_info);

%% functions
function A = largestNumericMatrix(S)
    names = fieldnames(S);
    bestName = '';
    bestNum  = -inf;
    for k = 1:numel(names)
        v = S.(names{k});
        if isnumeric(v) && ~isscalar(v) && ~isempty(v)
            if numel(v) > bestNum
                bestNum  = numel(v);
                bestName = names{k};
            end
        end
    end
    if bestNum < 0
        error('No numeric arrays found in the .mat file.');
    end
    A = S.(bestName);
end

% Vertical violin (posterior on y-axis, parameter index on x-axis)
% Vertical violin (posterior on y-axis, parameter index on x-axis)
function drawViolinV(x, x0, width, col, alphaFill, halfEye)
    x = x(:);
    x = x(x > 0);  
    if isempty(x)
        return;
    end

    
    logmin = log10(min(x));
    logmax = log10(max(x));
    dlog   = logmax - logmin;
    logmin = logmin - 0.05*dlog; 
    logmax = logmax + 0.05*dlog;

    yy = logspace(logmin, logmax, 256); 
    [f, yi] = ksdensity(x, yy);
    f = f / max(f);
    half = width * f;

    if halfEye
        % half-eye: only to the right of the center
        X = [x0*ones(size(half)), x0 + fliplr(half)];
        Y = [yi,                  fliplr(yi)];
    else
        % full violin
        X = [x0 - half, x0 + fliplr(half)];
        Y = [yi,        fliplr(yi)];
    end

    p = patch('XData',X,'YData',Y, ...
              'FaceColor',col,'EdgeColor','none','FaceAlpha',alphaFill);
    uistack(p,'bottom');
end


function plot_aoi_interval_compare_defau(build_sim, build_obs,build_sim_de, road_sim, road_obs, road_sim_de, aoi_info)
% Visualise comparisons between estimated and observed values across 25 AOIs
%  - x-axis: AOI (sorted by observed number of collapsed buildings, descending)
%  - y-axis: counts (collapsed buildings / closed road segments)
%  - For each AOI, plot:
%       * 2.5%–97.5% prediction interval (vertical line)
%       * median (marker)
%       * observed value (hollow diamond)

    % ------------ AOI labels (use IDs from aoi_info) -------------
    if nargin < 5 || isempty(aoi_info)
        NA = numel(build_obs);
        aoi_labels = arrayfun(@(i) sprintf('AOI %d', i), 1:NA, 'UniformOutput', false);
    else
        ids = aoi_info(:,1);
        aoi_labels = cellstr(num2str(ids));  % string ID for each AOI
    end

    % ------------ Standardise matrix dimensions: [NAOI x Nsim] ----------------
    NA = numel(build_obs);

    % building
    if size(build_sim,1) == NA
        S_build = build_sim;
    elseif size(build_sim,2) == NA
        S_build = build_sim.';
    else
        error('build_sim dimensions do not match build_obs');
    end

    % road
    if size(road_sim,1) == NA
        S_road = road_sim;
    elseif size(road_sim,2) == NA
        S_road = road_sim.';
    else
        error('road_sim dimensions do not match road_obs');
    end

    % building (default)
    if size(build_sim_de,1) == NA
        S_build_de = build_sim_de;
    elseif size(build_sim_de,2) == NA
        S_build_de = build_sim_de.';
    else
        error('build_sim dimensions do not match build_obs');
    end

    % road (default)
    if size(road_sim_de,1) == NA
        S_road_de = road_sim_de;
    elseif size(road_sim_de,2) == NA
        S_road_de = road_sim_de.';
    else
        error('road_sim dimensions do not match road_obs');
    end


    % ------------ Compute 2.5 / 50 / 97.5 percentiles -----------------
    % Each row corresponds to one AOI
    qB = prctile(S_build, [2.5 50 97.5], 2);  % [NA x 3]
    loB  = qB(:,1);
    medB = qB(:,2);
    hiB  = qB(:,3);

    qR = prctile(S_road, [2.5 50 97.5], 2);
    loR  = qR(:,1);
    medR = qR(:,2);
    hiR  = qR(:,3);

    qB_de = prctile(S_build_de, [2.5 50 97.5], 2);  % [NA x 3]
    loB_de  = qB_de(:,1);
    medB_de = qB_de(:,2);
    hiB_de  = qB_de(:,3);

    qR_de = prctile(S_road_de, [2.5 50 97.5], 2);
    loR_de  = qR_de(:,1);
    medR_de = qR_de(:,2);
    hiR_de  = qR_de(:,3);


    % ------------ Sort by observed collapsed buildings (descending) -----------------
    % [~, order] = sort(build_obs, 'descend');
    [~, order] = sortrows([build_obs,medB], 'descend');
    loB  = loB(order);  medB = medB(order);  hiB  = hiB(order);
    obsB = build_obs(order);
    loB_de  = loB_de(order);  medB_de = medB_de(order);  hiB_de  = hiB_de(order);

    % Reorder
    % [~, order] = sort(road_obs, 'descend');
    % [~, order] = sortrows([road_obs,medR], 'descend');
    loR  = loR(order);  medR = medR(order);  hiR  = hiR(order);
    loR_de  = loR_de(order);  medR_de = medR_de(order);  hiR_de  = hiR_de(order);
    obsR = road_obs(order);
    labels_sorted = aoi_labels(order);

    x = 1:NA;

    % Colours
    col_build = [0.20 0.45 0.85];
    col_road  = [0.90 0.50 0.10];
    col_ci    = [0.7 0.7 0.7];

    % ------------ Plot ----------------------------------------
    figure('Color','w','Position',[100 100 900 600]);

    % ----- Subplot 1: collapsed buildings -----
    subplot(1,2,1); hold on;

    % default
    % 95% interval: vertical lines
    for i = 1:NA
        plot([x(i) x(i)], [loB_de(i) hiB_de(i)], '-', ...
             'Color', [0.7,0.7,0.7], 'LineWidth', 3);
    end

    % Median: marker
    plot(x, medB_de, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.7,0.7,0.7], ...
         'MarkerEdgeColor', [0.7,0.7,0.7]);
         % 'MarkerFaceColor', [0.7,0.7,0.7], ...

    % calibrated
    % 95% interval: vertical lines
    for i = 1:NA
        plot([x(i) x(i)], [loB(i) hiB(i)], '-', ...
             'Color', col_build, 'LineWidth', 1.5);
    end

    % Median: marker
    plot(x, medB, 'x', 'MarkerSize', 9, ...
         'MarkerFaceColor', col_build, ...
         'MarkerEdgeColor', col_build);

    % Observed values: hollow diamonds
    plot(x, obsB, 'd', 'MarkerSize', 7, ...
         'MarkerFaceColor', 'none', ...
         'MarkerEdgeColor', col_road, ...
         'LineWidth', 1.2);

    xlim([0.5 NA+0.5]);
    ylabel('Collapsed buildings');
    xlabel('AOI');
    % title('AOI-level collapsed buildings: posterior vs observed');
    % set(gca,'XTick', x, 'XTickLabel', labels_sorted);
    grid off; box off;

    % legend({'95% posterior interval','Posterior median','Observed'}, ...
    %        'Location','northoutside','Orientation','horizontal');

    % ----- Subplot 2: closed road segments -----
    subplot(1,2,2); hold on;


    % default
    % 95% interval: vertical lines
    for i = 1:NA
        plot([x(i) x(i)], [loR_de(i) hiR_de(i)], '-', ...
             'Color', [0.7,0.7,0.7], 'LineWidth', 3);
    end

    % Median: marker
    plot(x, medR_de, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.7,0.7,0.7], ...
         'MarkerEdgeColor', [0.7,0.7,0.7]);
         % 'MarkerFaceColor', [0.7,0.7,0.7], ...
    
    % calibrated
    for i = 1:NA
        plot([x(i) x(i)], [loR(i) hiR(i)], '-', ...
             'Color', col_build, 'LineWidth', 1.5);
    end

    plot(x, medR, 'x', 'MarkerSize', 9, ...
         'MarkerFaceColor', col_build, ...
         'MarkerEdgeColor', col_build);

    plot(x, obsR, 'd', 'MarkerSize', 7, ...
         'MarkerFaceColor', 'none', ...
         'MarkerEdgeColor', col_road, ...
         'LineWidth', 1.2);

    xlim([0.5 NA+0.5]);
    xlabel('AOI');
    ylabel('Closed road segments');
    % title('AOI-level closed roads: posterior vs observed');
    % set(gca,'XTick', x, 'XTickLabel', labels_sorted);
    grid off; box off;

end
