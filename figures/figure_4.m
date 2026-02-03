%% 4a & 4b
clc
clear

result_dir = "figures\results";

city_ids = 1:15;
city_names = {'Ayeyarwady'	'Bago'	'Chin'	'Kachin'	'Kayah'	'Kayin'	'Magway'	'Mandalay'	'Mon'	'Naypyitaw'	'Rakhine'	'Sagaing'	'Shan'	'Tanintharyi'	'Yangon'};
cities_shp = shaperead(fullfile(result_dir,"MMR-1.shp"));
gids_num = cellfun(@(x) sscanf(x, 'MMR.%d_1'), {cities_shp.GID_1});

for ci = 1:15
    GID_code = city_ids(ci);
    city_shps(ci) = cities_shp(gids_num==GID_code);
end

load(fullfile(result_dir,"all_cell_risk_results.mat"), ...
    "all_cell_loss_post","all_cell_fatality_post", ...
    "cell_informs_uni","all_cell_loss_base","all_cell_fatality_base")

load(fullfile(result_dir,"gid1_risk.mat"),"GID1_sum");
GID1_post_AAL = GID1_sum.Build_Loss_post;
GID1_post_AAF = GID1_sum.Fatality_post;

% 4a
out_AAF = plot_fig3_grid_plus_regions( ...
    cell_informs_uni, ...
    all_cell_fatality_post, ...
    GID1_post_AAF, ...
    city_shps, ...
    city_names, ...
    'Units','F', ...
    'CellPercentileCuts',[1 5 10 25], ...              % your custom bands
    'RegionValueEdges',[0 1 2 5 10 Inf], ...  % exact edges (5 classes)
    'RegionSizeRange',[24 160], 'RegionScale',0.9);


% 4b
out_AAL = plot_fig3_grid_plus_regions( ...
    cell_informs_uni, ...
    all_cell_loss_post, ...
    GID1_post_AAL, ...
    city_shps, ...
    city_names, ...
    'Units','USD', ...
    'CellPercentileCuts',[1 5 10 25], ...              % your custom bands
    'RegionValueEdges',[0 1.0e+05 1.0e+06 1.0e+07 5.0e+07 Inf], ...  % exact edges (5 classes)
    'RegionSizeRange',[24 160], 'RegionScale',0.9);

%% 4d & 4f
clear
% 1) Load required data (GEM shapefile) and assemble variables used in the figure
result_dir = "figures\results";

gem_shp_path = fullfile(result_dir,"mmr_gem_risk_map_with_cali_default.shp");
gem_shp = shaperead(gem_shp_path);

gem_data.lon_lat     = [[gem_shp.lon]', [gem_shp.lat]'];
gem_data.losses      = [gem_shp.losses]';
gem_data.fatalities  = [gem_shp.fatalities]';
gem_data.buildings   = [gem_shp.buildings]';
gem_data.area        = [gem_shp.area]';

gem_data.cali_losses    = [gem_shp.loss_sum]';
gem_data.cali_fata      = [gem_shp.fata_sum]';
gem_data.cali_deflosses = [gem_shp.deflos__su]';
gem_data.cali_deffata   = [gem_shp.deffat_sum]';

% 2) Map each point to a region (GID3) and filter out points without a valid region
gid3BoundaryFile = fullfile(result_dir,"MMR-3.shp");
[gid, ~] = map_building_gid(gem_data.lon_lat, gid3BoundaryFile);

has_gid = gid(:,1) > 0;

% Filter all fields so lengths remain consistent
gem_data.losses        = gem_data.losses(has_gid,:);
gem_data.fatalities    = gem_data.fatalities(has_gid,:);
gem_data.buildings     = gem_data.buildings(has_gid,:);
gem_data.area          = gem_data.area(has_gid,:);

gem_data.cali_losses    = gem_data.cali_losses(has_gid,:);
gem_data.cali_fata      = gem_data.cali_fata(has_gid,:);
gem_data.cali_deflosses = gem_data.cali_deflosses(has_gid,:);
gem_data.cali_deffata   = gem_data.cali_deffata(has_gid,:);

gid = gid(has_gid,:);

% 3) Figure: regional share of top 10% hotspots, GEM vs calibrated (fatalities + losses)
city_names = {'Ayeyarwady' 'Bago' 'Chin' 'Kachin' 'Kayah' 'Kayin' 'Magway' ...
              'Mandalay' 'Mon' 'Naypyitaw' 'Rakhine' 'Sagaing' 'Shan' ...
              'Tanintharyi' 'Yangon'};

region_id = double(gid(:,1));
nRegion   = numel(city_names);
N         = numel(gid(:,1));

topFrac = 0.1;
K = round(topFrac * N);

figure;

% (d) Fatalities: GEM vs calibrated
loss_base = gem_data.fatalities;  loss_base(isnan(loss_base)) = 0;   % GEM
loss_post = gem_data.cali_fata;     loss_post(isnan(loss_post)) = 0;   % calibrated pipeline

[~, idx_sorted_base] = sort(loss_base, 'descend');
idx_top_base = idx_sorted_base(1:K);

[~, idx_sorted_post] = sort(loss_post, 'descend');
idx_top_post = idx_sorted_post(1:K);

count_base = accumarray(region_id(idx_top_base), 1, [nRegion 1]);
count_post = accumarray(region_id(idx_top_post), 1, [nRegion 1]);

% Sort regions once (and reuse for the second panel)
[~, sort_idx] = sort(count_base, 'ascend');

subplot(1,2,1); hold on
[share_base_10_aaf, share_post_10_aaf] = plot_region_share_dumbbell( ...
    city_names, count_base, count_post, sort_idx, ...
    'Share of top-1% loss hotspots (%)', ...
    'Regional share of top-1% loss hotspots: GEM vs calibrated');

% (f) Losses: GEM vs calibrated
loss_base = gem_data.losses;  loss_base(isnan(loss_base)) = 0; % GEM
loss_post = gem_data.cali_losses;     loss_post(isnan(loss_post)) = 0; % calibrated pipeline

[~, idx_sorted_base] = sort(loss_base, 'descend');
idx_top_base = idx_sorted_base(1:K);

[~, idx_sorted_post] = sort(loss_post, 'descend');
idx_top_post = idx_sorted_post(1:K);


count_base = accumarray(region_id(idx_top_base), 1, [nRegion 1]);
count_post = accumarray(region_id(idx_top_post), 1, [nRegion 1]);

subplot(1,2,2); hold on
[share_base_10_aal, share_post_10_aal] = plot_region_share_dumbbell( ...
    city_names, count_base, count_post, sort_idx, ...
    'Share of top-1% loss hotspots (%)', ...
    'Regional share of top-1% loss hotspots: GEM vs calibrated');


%% functions

%
function out = plot_fig3_grid_plus_regions(cell_informs_uni, cell_values, region_sims, city_shps, city_names, varargin)
% PLOT_FIG3_GRID_PLUS_REGIONS (discrete classes with user-specified cuts)
% Layer order: base ADM1 fill (light gray) -> grid (discrete) -> ADM1 boundary -> region circles (discrete)
%
% Inputs
%   cell_informs_uni : [N x 3] = [cell_id, lon, lat]
%   cell_values      : [N x 1] median per cell (>=0)
%   region_sims      : [R x M] or [M x R] posterior draws (only median used for classes)
%   city_shps        : 1xR shaperead structs (X=lon, Y=lat)
%   city_names       : 1xR cellstr
%
% Name-Value options
%   'Basemap'            : 'grayland'
%   'BaseColor'          : [0.95 0.95 0.95]
%   'DrawBaseFill'       : true
%   'FigTitle'           : 'Fig. 3 — Estimated impact (2025 event)'
%   'Units'              : ''       % e.g., 'fatalities' | 'USD'
%   'CellMarkerSize'     : 6
%   'AdminEdgeColor'     : [0 0 0]
%   'AdminLineWidth'     : 0.6
%   'ShowLabels'         : true
%   'RegionSizeRange'    : [30 220] % circle area (points^2)
%   'RegionScale'        : 1.0
%   'CellPercentileCuts' : [1 5 10 20]  % top-% cuts (ascending), producing 5 classes
%   'RegionValueEdges'   : []            % exact numeric edges for region medians; e.g., [0 500 2000 5000 Inf]
%   'SaveAs'             : ''
%
% Output
%   out.cell_percentile_cuts : vector of supplied percentile cuts
%   out.cell_value_thresholds: corresponding value thresholds (ascending, complements of cuts)
%   out.region_edges         : final region edges used (ascending, length=K+1)
%   out.region_median        : [R x 1] medians
%   out.centroids            : [R x 2] [lon, lat]

% ---------------- Parse options ----------------
p = inputParser;
addParameter(p,'Basemap','grayland');
addParameter(p,'BaseColor',[0.95 0.95 0.95]);
addParameter(p,'DrawBaseFill',true,@islogical);
addParameter(p,'FigTitle','');
addParameter(p,'Units','');
addParameter(p,'CellMarkerSize',6);
addParameter(p,'AdminEdgeColor',[0 0 0]);
addParameter(p,'AdminLineWidth',0.6,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'ShowLabels',true,@islogical);
addParameter(p,'RegionSizeRange',[30 220]);
addParameter(p,'RegionScale',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'CellPercentileCuts',[1 5 10 20],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'RegionValueEdges',[],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'SaveAs','');
addParameter(p,'LegendTokenPad',1.6,@(x)isnumeric(x)&&isscalar(x)&&x>0);

parse(p,varargin{:});
opt = p.Results;

% ---------------- Data prep ----------------
lon = cell_informs_uni(:,2);
lat = cell_informs_uni(:,3);
v   = max(cell_values(:),0);

% Region medians
if size(region_sims,1) == numel(city_names), RS = region_sims; else, RS = region_sims.'; end
medR = median(RS,2,'omitnan');

R = numel(city_shps);
cent = nan(R,2);
for i = 1:R
    xi = city_shps(i).X(:); yi = city_shps(i).Y(:);
    m  = ~isnan(xi) & ~isnan(yi);
    cent(i,:) = [mean(xi(m)), mean(yi(m))];   % [lon, lat]
end

% ---------------- Figure & basemap ----------------
f = figure('Color','w','Position',[80 80 300 840]);
gx = geoaxes('Parent',f); hold(gx,'on');
try geobasemap(gx, string(opt.Basemap)); catch, geobasemap(gx,'grayland'); end
try grid(gx,'on'); catch, grid on; end

% ---------------- BASE FILL (light gray) ----------------
if opt.DrawBaseFill
    try
        polys = geopolyshape.empty;
        for i=1:R, polys(i) = geopolyshape(city_shps(i).Y, city_shps(i).X); end
        hp = geoplot(gx, polys, 'FaceColor', opt.BaseColor, 'EdgeColor', 'none');
        uistack(hp,'bottom');
    catch
        % fallback silently
    end
end

% ---------------- CELLS: discrete classes by percentile cuts ----------------
% Normalize & validate cuts: ascending, within (0,100)
pc = sort(unique(opt.CellPercentileCuts(:)'));
pc = pc(pc>0 & pc<100);
if isempty(pc), pc = [1 5 10 20]; end
% Complement percentiles → quantiles (ascending): e.g., [80 90 95 99]
qpr = sort(100 - pc);
% Value thresholds at those quantiles
qVals = prctile(v, qpr);   % ascending thresholds
% Build class masks: other, then bands up to top
classMasks = cell(numel(qVals)+1,1);
classMasks{1} = v < qVals(1);                    % other (< qVal1)
for j = 1:numel(qVals)-1
    classMasks{j+1} = v >= qVals(j) & v < qVals(j+1);
end
classMasks{end} = v >= qVals(end);                % top (>= last)
nClasses = numel(classMasks);

% Colors: very light gray for "other", then gradient from pale orange to dark red
cols = cellClassColors(nClasses);

% Draw cells from low to high class to keep hotspots on top
msz = opt.CellMarkerSize; edgeC = [0 0 0]; lw = 0.05;
for j = 1:nClasses
    idx = classMasks{j};
    if any(idx)
        geoscatter(gx, lat(idx), lon(idx), msz, cols(j,:), 's','filled','MarkerEdgeColor',edgeC,'LineWidth',lw);
    end
end

% View
geolimits(gx, [min(lat)-0.2 max(lat)+0.2], [min(lon)-0.2 max(lon)+0.2]);

% ---------------- ADMIN boundaries (above cells, below circles) ----------------
for i=1:R
    geoplot(gx, city_shps(i).Y, city_shps(i).X, 'Color', opt.AdminEdgeColor, 'LineWidth', opt.AdminLineWidth);
end

% ---------------- REGIONS: discrete classes by exact value edges ----------------
if ~isempty(opt.RegionValueEdges)
    edges = sort(opt.RegionValueEdges(:)');
    % Ensure coverage of data range
    if edges(1) > -Inf, edges = [min(min(medR),edges(1)) edges]; end
    if edges(end) <  Inf, edges = [edges max(max(medR),edges(end))]; end
else
    % Default to quantile-based edges (5 classes)
    Kdef = 5;
    edges = prctile(medR, linspace(0,100,Kdef+1));
    edges(1) = min(edges(1), min(medR)); edges(end) = max(edges(end), max(medR));
end
% Clean duplicates to avoid empty/ill-defined bins
edges = unique(edges,'stable');
if numel(edges) < 2, edges = [min(medR) max(medR)]; end
K = numel(edges)-1;

Amin = opt.RegionSizeRange(1); Amax = opt.RegionSizeRange(2);
size_per_class = linspace(Amin, Amax, K) * opt.RegionScale;
circFace = [0.98 0.80 0.25];  circEdge = [0.35 0.35 0.35];

% Draw circles per class (increasing size with class index)
for k = 1:K
    if k < K
        inK = (medR >= edges(k)) & (medR < edges(k+1));
        labK = sprintf('%s – %s', numfmt(edges(k)), numfmt(edges(k+1)));
    else
        inK = (medR >= edges(k)) & (medR <= edges(k+1));
        labK = sprintf('\x2265 %s', numfmt(edges(k)));
    end
    if any(inK)
        geoscatter(gx, cent(inK,2), cent(inK,1), size_per_class(k), ...
            'o','MarkerFaceColor',circFace,'MarkerEdgeColor',circEdge,'LineWidth',0.8, ...
            'DisplayName', ['Region median: ' labK]);
    end
end

% Labels (optional, only median)
if opt.ShowLabels
    for i = 1:R
        lbl = sprintf('%s — %s', string(city_names{i}), numfmt(medR(i)));
        text(gx, cent(i,1), cent(i,2), lbl, ...
            'FontSize',8,'Color',[0 0 0],'BackgroundColor',[1 1 1 0.85], ...
            'Margin',3,'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

% ---------------- Combined legend (cells + regions) ----------------
unitStr = opt.Units; if ~isempty(unitStr), unitStr = [' ' unitStr]; end

% Cells legend proxies & labels (from top to low)
% Build human-readable bands using complements of qpr
topBands = 100 - qpr;           % e.g., [20 10 5 1] for qpr=[80 90 95 99]
nCuts = numel(qVals);
h_cell = gobjects(nClasses,1); t_cell = strings(nClasses,1);

% Top band
h_cell(1) = plot(nan,nan,'s','MarkerSize',9,'MarkerFaceColor',cols(end,:),'MarkerEdgeColor','k');
t_cell(1) = sprintf('Top %g%% (\x2265 %s%s)', topBands(end), numfmt(qVals(end)), unitStr);
% Middle bands (descending importance)
for j = nCuts-1:-1:1
    ci = nCuts - j + 1;                   % legend row index
    h_cell(ci) = plot(nan,nan,'s','MarkerSize',9,'MarkerFaceColor',cols(j+1,:),'MarkerEdgeColor','k');
    loV = qVals(j); hiV = qVals(j+1);
    loP = topBands(j+1); hiP = topBands(j);
    t_cell(ci) = sprintf('%g–%g%% ([%s–%s)%s)', loP, hiP, numfmt(loV), numfmt(hiV), unitStr);
end
% Other
h_cell(end) = plot(nan,nan,'s','MarkerSize',9,'MarkerFaceColor',cols(1,:),'MarkerEdgeColor','k');
t_cell(end) = sprintf('Other (< %s%s)', numfmt(qVals(1)), unitStr);

% Region legend proxies (increasing size)
h_reg = gobjects(K,1); t_reg = strings(K,1);
for k = 1:K
    h_reg(k) = plot(nan,nan,'o','MarkerSize',sqrt(size_per_class(k)), ...
        'MarkerFaceColor',circFace,'MarkerEdgeColor',circEdge,'LineWidth',0.8);
    if k < K
        t_reg(k) = sprintf('%s – %s%s', numfmt(edges(k)), numfmt(edges(k+1)), unitStr);
    else
        t_reg(k) = sprintf('\x2265 %s%s', numfmt(edges(k)), unitStr);
    end
end

% Assemble legend: cells (top→low) then regions
leg = legend(gx, [h_cell; h_reg], [t_cell; t_reg], 'Location','southoutside','Orientation','vertical');
% leg.ItemTokenSize = [12 12];

msz_max = 0;
for k = 1:numel(h_reg)
    if isgraphics(h_reg(k))
        msz_max = max(msz_max, get(h_reg(k),'MarkerSize'));   % MarkerSize is in points
    end
end
msz_max = max(msz_max, 9);  % at least as big as the cell squares

pad = opt.LegendTokenPad;   % e.g., 1.6 -> some headroom
try
    % ItemTokenSize = [token width, token height] in points
    leg.ItemTokenSize = [max(leg.ItemTokenSize(1), ceil(msz_max*pad)) ...
                         max(leg.ItemTokenSize(2), ceil(msz_max*pad))];
catch
    % Older MATLAB (no ItemTokenSize): shrink proxy markers to fit
    shrink = 0.75;  % gentle fallback
    for k = 1:numel(h_reg)
        if isgraphics(h_reg(k))
            set(h_reg(k),'MarkerSize',get(h_reg(k),'MarkerSize')*shrink);
        end
    end
end

title(gx, string(opt.FigTitle));

% Save
if ~isempty(opt.SaveAs)
    fn = char(opt.SaveAs);
    exportgraphics(f, fn, 'Resolution', 300);
    try
        [p,n,~]=fileparts(fn);
        exportgraphics(f, fullfile(p,[n '.pdf']), 'ContentType','vector');
    end
end

% Outputs
out = struct();
out.cell_percentile_cuts  = pc;
out.cell_value_thresholds = qVals;
out.region_edges          = edges;
out.region_median         = medR;
out.centroids             = cent;
end

% ===== Helpers =====
function cols = cellClassColors(n)
% First class = "other" very light gray; remaining form a gradient pale→dark
if n<=1, cols = [0.95 0.95 0.95]; return; end
c_other = [0.95 0.95 0.95];
c_lo    = [0.99 0.85 0.60];   % pale orange
c_hi    = [0.75 0.00 0.00];   % dark red
cols = zeros(n,3);
cols(1,:) = c_other;
for k = 2:n
    t = (k-2)/max(1,(n-2));   % 0..1 across non-"other" bands
    cols(k,:) = (1-t)*c_lo + t*c_hi;
end
end

function s = numfmt(x)
% NUMFMT  Magnitude-aware number formatter.
    if isnan(x) || isinf(x)
        s = sprintf('%g', x);
        return;
    end
    if x == 0
        s = '0';
        return;
    end

    ax = abs(x);

    if ax >= 1000
        fmt = java.text.DecimalFormat('###,##0');
        s = char(fmt.format(java.lang.Double(x)));
        return;
    end

    if ax >= 100
        s = sprintf('%.0f', x);
    elseif ax >= 10
        s = sprintf('%.1f', x);
    elseif ax >= 1
        s = sprintf('%.2f', x);
    elseif ax >= 0.01
        s = sprintf('%.3f', x);
    elseif ax >= 0.001
        s = sprintf('%.4f', x);
    else
        s = sprintf('%.2e', x);
    end
end

function [share_base,share_post]=plot_region_share_dumbbell(city_names, share_base, share_post, sort_idx,...
    yLabelStr, titleStr)
%PLOT_REGION_SHARE_DUMBBELL Vertical dumbbell for regional hotspot shares.
%
%   plot_region_share_dumbbell(city_names, share_base, share_post)
%   compares the regional shares of hotspots between a default and a
%   calibrated model. Inputs share_base and share_post can be counts or
%   percentages; they will be normalised to sum to 100%.
%
%   plot_region_share_dumbbell(..., yLabelStr, titleStr) lets you specify
%   the y-axis label and figure title.
%
%   city_names : 1xR or Rx1 cell array of region names
%   share_base : Rx1 numeric vector (default model)
%   share_post : Rx1 numeric vector (calibrated model)

    % --- input handling ---------------------------------------------------
    if nargin < 4 || isempty(yLabelStr)
        yLabelStr = 'Share of top-1% loss hotspots (%)';
    end
    if nargin < 5 || isempty(titleStr)
        titleStr  = 'Regional share of top-1% hotspots: default vs calibrated';
    end

    % ensure column vectors
    share_base = share_base(:);
    share_post = share_post(:);

    nRegion = numel(share_base);
    if numel(share_post) ~= nRegion
        error('share_base and share_post must have the same length.');
    end
    if numel(city_names) ~= nRegion
        error('city_names length must match share vectors.');
    end

    % normalise to percentages (each model sums to 100)
    share_base = 100 * share_base ./ sum(share_base);
    share_post = 100 * share_post ./ sum(share_post);

    % sort regions by calibrated share (descending)
    % [share_post_sorted, sort_idx] = sort(share_post, 'descend');
    % share_base_sorted = share_base(sort_idx);
    share_base_sorted = share_base(sort_idx);
    share_post_sorted = share_post(sort_idx);

    city_names_sorted = city_names(sort_idx);

    % --- plotting ---------------------------------------------------------
    x = 1:nRegion;

    % figure; hold on;

    for r = 1:nRegion
        base_r = share_base_sorted(r);
        post_r = share_post_sorted(r);

        % colour: red if calibrated share increased, blue if decreased
        if post_r > base_r
            col = [0.85 0.33 0.10];   % red
        elseif post_r < base_r
            col = [0.00 0.45 0.74];   % blue
        else
            col = [0.5 0.5 0.5];      % grey (no change)
        end

        % % vertical line: default -> calibrated
        % plot([x(r) x(r)], [base_r post_r], '-', ...
        %     'Color', col, 'LineWidth', 1);
        % 
        % % small hollow circle: default
        % plot(x(r), base_r, 'o', ...
        %     'MarkerSize', 4, ...
        %     'MarkerEdgeColor', col, ...
        %     'MarkerFaceColor', 'w');
        % 
        % % large filled circle: calibrated
        % plot(x(r), post_r, 'o', ...
        %     'MarkerSize', 7, ...
        %     'MarkerEdgeColor', col, ...
        %     'MarkerFaceColor', col);

        % vertical line: default -> calibrated
        plot([base_r post_r], [x(r) x(r)], '-', ...
            'Color', col, 'LineWidth', 1);

        % small hollow circle: default
        plot(base_r, x(r), 'o', ...
            'MarkerSize', 4, ...
            'MarkerEdgeColor', col, ...
            'MarkerFaceColor', 'w');

        % large filled circle: calibrated
        plot(post_r, x(r), 'o', ...
            'MarkerSize', 7, ...
            'MarkerEdgeColor', col, ...
            'MarkerFaceColor', col);

    end

    % xlim([0.5 nRegion+0.5]);
    % ylim([0, max([share_base_sorted; share_post_sorted]) * 1.1]);
    % 
    % set(gca, 'XTick', x, 'XTickLabel', city_names_sorted, ...
    %          'XTickLabelRotation', 45, ...
    %          'TickDir', 'out');

    ylim([0.5 nRegion+0.5]);
    xlim([0, max([share_base_sorted; share_post_sorted]) * 1.1]);

    set(gca, 'YTick', x, 'YTickLabel', city_names_sorted, ...
             'TickDir', 'out');

    % ylabel(yLabelStr);
    grid off;
    box off;

    % legend: small = default, large = calibrated
    % h_def = plot(nan, nan, 'o', 'MarkerSize', 4, ...
    %     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
    % h_cal = plot(nan, nan, 'o', 'MarkerSize', 7, ...
    %     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    % legend([h_def, h_cal], {'Default model', 'Calibrated model'}, ...
    %        'Location', 'northoutside', 'Orientation', 'horizontal');
    % 
    % title(titleStr);
end
