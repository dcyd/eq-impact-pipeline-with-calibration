%% Figs. 3a and 3b
clc
clear
result_dir = "figures\results";

load(fullfile(result_dir,"cell_impacts.mat"), ...
    "all_cell_fatality_post","all_cell_loss_post","cell_informs_uni");

load(fullfile(result_dir,"gid1_results_2025.mat"), ...
    "gid1_Fatality_post", ...
    "loss_post_gid1")

city_ids = 1:15;
city_names = {'Ayeyarwady'	'Bago'	'Chin'	'Kachin'	'Kayah'	'Kayin'	'Magway'	'Mandalay'	'Mon'	'Naypyitaw'	'Rakhine'	'Sagaing'	'Shan'	'Tanintharyi'	'Yangon'};
cities_shp = shaperead(fullfile(result_dir,"MMR-1.shp"));
gids_num = cellfun(@(x) sscanf(x, 'MMR.%d_1'), {cities_shp.GID_1});

for ci = 1:15
    GID_code = city_ids(ci);
    city_shps(ci) = cities_shp(gids_num==GID_code);
end

% fig. 3a

plot_fig3_grid_plus_regions( ...
    cell_informs_uni, ...
    all_cell_fatality_post, ...
    gid1_Fatality_post, ...
    city_shps, ...
    city_names, ...
    'Units','F', ...
    'CellPercentileCuts',[1 5 10 25], ...              % your custom bands
    'RegionValueEdges',[0 200 500 1000 2500 Inf], ...  % exact edges (5 classes)
    'RegionSizeRange',[24 160], 'RegionScale',0.9);

% fig. 3b
plot_fig3_grid_plus_regions( ...
    cell_informs_uni, ...
    all_cell_loss_post, ...
    loss_post_gid1, ...
    city_shps, ...
    city_names, ...
    'Units','USD', ...
    'CellPercentileCuts',[1 5 10 25], ...              % your custom bands
    'RegionValueEdges',[0 1.0e+07 1.0e+08 5.0e+08 1.0e+09 Inf], ...  % exact edges (5 classes)
    'RegionSizeRange',[24 160], 'RegionScale',0.9);

%% Figs. 3c-cf
clear

result_dir = "figures\results";
load(fullfile(result_dir,"gid1_results_2025.mat"),"gid1_Fatality_post","loss_post_gid1");

gid0_building_loss_post_2025 = sum(loss_post_gid1);
gid0_fatality_post_2025 = sum(gid1_Fatality_post);

load(fullfile(result_dir,"gid1_results_2012.mat"),"gid1_Fatality_post");
gid0_fatality_post_2012 = sum(gid1_Fatality_post);

load(fullfile(result_dir,"gid1_results_2016.mat"),"gid1_Fatality_post");
gid0_fatality_post_2016 = sum(gid1_Fatality_post);

gid1_building_loss_obs = [92 918 0 0 12 20 306 3686 20 631 0 1697 154 0 69] * 1e6;
gid0_building_loss_obs_2025 = sum(gid1_building_loss_obs);
gid0_fatality_obs_2025 = 3757;

gid0_fatality_obs_2012 = 38;
gid0_fatality_obs_2016 = 4;

plot_gid0_loss_posterior(gid0_building_loss_post_2025/1000000, gid0_building_loss_obs_2025/1000000)
xlabel('L(2025), million');

plot_gid0_loss_posterior(gid0_fatality_post_2025, gid0_fatality_obs_2025)
xlabel('F(2025)');

plot_gid0_loss_posterior(gid0_fatality_post_2016, gid0_fatality_obs_2016)
xlabel('F(2016)');

plot_gid0_loss_posterior(gid0_fatality_post_2012, gid0_fatality_obs_2012)
xlabel('F(2012)');


%% functions

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


function S = plot_gid0_loss_posterior(post_samples, obs_value, varargin)
% Non-negative-domain version (axes swapped): X = value (>=0), Y = density (hidden)
% Shows a vertical histogram + KDE + vertical interval bands for 80/90/95%
% + vertical lines for median / observed value.
%
% Name-Value:
%   'Title'  : title
%   'Bins'   : number of histogram bins (default: automatic)
%   'LogY'   : log scale on the "value axis" (now the X axis; default false)
%   'SaveAs' : output filename for exporting

% ---------- Parameters ----------
p = inputParser;
addParameter(p,'Title','GID0 building loss posterior',@(x)ischar(x)||isstring(x));
addParameter(p,'Bins',[],@(x) isempty(x) || (isscalar(x)&&x>0));
addParameter(p,'LogY',false,@islogical);   % keep the name, but it controls the X axis
addParameter(p,'SaveAs','',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
ttl    = p.Results.Title;
nbins  = p.Results.Bins;
useLog = p.Results.LogY;   % now means: use log scale on X=loss axis
saveAs = p.Results.SaveAs;

% ---------- Data and quantiles (clamped to non-negative) ----------
x   = max(post_samples(:), 0);     % non-negative
obs = max(double(obs_value), 0);   % non-negative

S.med     = median(x);
S.pi80_lo = max(0, prctile(x,10));  S.pi80_hi = prctile(x,90);
S.pi90_lo = max(0, prctile(x,5));   S.pi90_hi = prctile(x,95);
S.pi95_lo = max(0, prctile(x,2.5)); S.pi95_hi = prctile(x,97.5);
S.obs     = obs;
S.in80    = (obs>=S.pi80_lo)&&(obs<=S.pi80_hi);
S.in90    = (obs>=S.pi90_lo)&&(obs<=S.pi90_hi);
S.in95    = (obs>=S.pi95_lo)&&(obs<=S.pi95_hi);

% ---------- Number of histogram bins ----------
if isempty(nbins)
    n = numel(x);
    if n>1
        h = 2*iqr(x)/(n^(1/3));
        if h>0
            nbins = max(10, ceil((max(x)-0)/h)); % start from 0
        else
            nbins = 30;
        end
    else
        nbins = 30;
    end
end

% ---------- Plot ----------
figure('Color','w','Position',[120 100 860 540]);
ax = axes; hold(ax,'on');

% Vertical histogram (X=value, starting from 0)
if ~useLog
    histogram(ax, x, nbins, 'Normalization','pdf', ...
        'FaceColor',[0.75 0.83 0.98], 'EdgeColor','none', ...
        'DisplayName','Posterior histogram', ...
        'BinLimits',[0 max(x)]);
else
    % Log X: cannot include 0, so use the smallest positive value as the lower bound
    pos = x(x>0);
    if isempty(pos), pos = 1; end
    xmin = min(pos);
    xmax = max(x);
    edges = logspace(log10(xmin), log10(max(xmax, xmin)), nbins+1);
    histogram(ax, x, edges, 'Normalization','pdf', ...
        'FaceColor',[0.75 0.83 0.98], 'EdgeColor','none', ...
        'DisplayName','Posterior histogram');
    set(ax,'XScale','log');
end

% KDE: positive support + boundary reflection correction
% (plot xi on X, f on Y)
try
    [f, xi] = ksdensity(x, 'Support','positive','BoundaryCorrection','reflection');
catch
    % Compatibility fallback for older MATLAB: estimate in log1p space (shape only)
    z  = log1p(x);
    zg = linspace(min(z), max(z), 512);
    [fz, zz] = ksdensity(z, zg);
    xi = exp(zz) - 1;     % >=0
    f  = fz / trapz(xi, fz);
end
plot(ax, xi, f, 'b-', 'LineWidth', 2.0, 'DisplayName','KDE');

% Current Y height (density axis) for filling interval bands
yhi = ylim(ax);
set(ax,'YLim',yhi);

% When using log X, PI lower bounds cannot be 0
x_eps = 0;
if useLog
    pos = x(x>0);
    if isempty(pos), pos = 1; end
    x_eps = max(min(pos)*0.9, 1e-12);
end

% Draw 95/90/80% vertical interval bands (light to dark), stacked at the bottom
p95 = plotPIBandV(ax, max(S.pi95_lo, x_eps), max(S.pi95_hi, x_eps), yhi(2), [0.80 0.80 0.80], 0.25);
p90 = plotPIBandV(ax, max(S.pi90_lo, x_eps), max(S.pi90_hi, x_eps), yhi(2), [0.60 0.60 0.60], 0.18);
p80 = plotPIBandV(ax, max(S.pi80_lo, x_eps), max(S.pi80_hi, x_eps), yhi(2), [0.30 0.30 0.30], 0.12);
uistack([p95 p90 p80],'bottom');

% Median and observed value: vertical lines (avoid drawing at 0 in log scale)
med_line = S.med;
obs_line = obs;
if useLog
    med_line = max(med_line, x_eps);
    obs_line = max(obs_line, x_eps);
end
xline(ax, med_line, 'k-', 'LineWidth',1.8, 'DisplayName','Median');
xline(ax, obs_line,  'r-', 'LineWidth',1.8, 'DisplayName','Observed');

% Axis styling: hide the density axis (Y axis)
try
    ax.YAxis.Visible = 'off';
catch
    set(ax,'YTick',[],'YColor',[1 1 1]);
end
box off; grid off; set(gca,'GridLineStyle',':');

end

% ===== Vertical interval band: fill the entire density axis (Y) for X in [lo, hi] =====
function h = plotPIBandV(ax, lo, hi, yhi, gray, alpha)
    xl = xlim(ax);
    lo = max(lo, 0);  % non-negative
    hi = max(hi, lo); % avoid hi < lo
    yl = ylim(ax);
    X = [lo hi hi lo];
    Y = [yl(1) yl(1) yhi yhi];
    h = patch('XData',X,'YData',Y,'FaceColor',gray, ...
              'EdgeColor','none','FaceAlpha',alpha,'Parent',ax);
    xlim(ax, xl);
end
