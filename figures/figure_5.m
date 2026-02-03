%% Figs. 5a and 5b
clc
clear

city_names = {'Ayeyarwady'	'Bago'	'Chin'	'Kachin'	'Kayah'	'Kayin'	'Magway'	'Mandalay'	'Mon'	'Naypyitaw'	'Rakhine'	'Sagaing'	'Shan'	'Tanintharyi'	'Yangon'};

result_dir = "figures\results";
load(fullfile(result_dir,"gid1_accessibiity_results_withbase.mat"),...
                "gid1_unacc_injury_60_post", ...
                "gid1_unacc_injury_60_init",...
                "gid1_unacc_injury_60_base", ...
                "gid1_unacc_injury_60_base_init");

gid1_unacc_injury_60_post = squeeze(gid1_unacc_injury_60_post(1,:,:));
gid1_unacc_injury_60_base = squeeze(gid1_unacc_injury_60_base(1,:,:));
gid1_unacc_injury_60_base_init = squeeze(gid1_unacc_injury_60_base_init(1,:,:));

gid0_unacc_injury_60_post_2025 =  sum(gid1_unacc_injury_60_post,1);
gid0_unacc_injury_60_init_2025 =  sum(gid1_unacc_injury_60_init,1);
gid0_unacc_injury_60_base_2025 =  sum(gid1_unacc_injury_60_base,1);
gid0_unacc_injury_60_base_init_2025 =  sum(gid1_unacc_injury_60_base_init,1);


gid0_unacc_injury_60 = [gid0_unacc_injury_60_post_2025;
            gid0_unacc_injury_60_init_2025;
            gid0_unacc_injury_60_base_2025;
            gid0_unacc_injury_60_base_init_2025;];
%
% OUT_impact = plot_four_pi_bars(gid0_unacc_injury_60, 'injuries with t_h>1h');
% xlim([0.5,4.5]);


%% Fig.5a

[~, ord] = sort(median(gid1_unacc_injury_60_post,  2,'omitnan'),'descend');

fig = figure('Color','w');
ax  = axes(fig); hold(ax,'on');

% ---------- medians per region
median_unacc  = [median(gid1_unacc_injury_60_post,  2,'omitnan'),...
                median(gid1_unacc_injury_60_init,  2,'omitnan'),...
                median(gid1_unacc_injury_60_base,  2,'omitnan'),...
                median(gid1_unacc_injury_60_base_init,  2,'omitnan'),];
x_lables = {'Impaired (calibrated)','Intact (calibrated)','Impaired (default)','Intact (default)'};
% x_lables = {'Impaired','Intact','Impaired','Intact'};

OUT = plot_unacc_injury_stacked_nature(median_unacc, ...
    city_names, ord, 'ShowTitle', false,...
    'labels',x_lables);

ax.FontName   = 'Helvetica';
ax.FontSize   = 8;
ax.LineWidth  = 0.75;
ax.TickDir    = 'out';
xlim([0.5,4.5])

%%
load(fullfile(result_dir,"risk_gid0_gid1_acc_results.mat"), ...
    "GID1_sum");

gid1_unacc_injury_60_post_risk = squeeze(GID1_sum.UnaccInj60_post(1,:,:));
gid1_unacc_injury_60_init_risk = GID1_sum.UnaccInj60_init;
gid1_unacc_injury_60_base_risk = squeeze(GID1_sum.UnaccInj60_base(1,:,:));
gid1_unacc_injury_60_base_init_risk = GID1_sum.UnaccInj60_base_init;

gid0_unacc_injury_60_post_risk =  sum(gid1_unacc_injury_60_post_risk,1);
gid0_unacc_injury_60_init_risk =  sum(gid1_unacc_injury_60_init_risk,1);
gid0_unacc_injury_60_base_risk =  sum(gid1_unacc_injury_60_base_risk,1);
gid0_unacc_injury_60_base_init_risk =  sum(gid1_unacc_injury_60_base_init_risk,1);

gid0_unacc_injury_60_risk = [gid0_unacc_injury_60_post_risk;
            gid0_unacc_injury_60_init_risk;
            gid0_unacc_injury_60_base_risk;
            gid0_unacc_injury_60_base_init_risk;];

% OUT_risk = plot_four_pi_bars(gid0_unacc_injury_60_risk, 'AA injuries with t_h>1h');
% xlim([0.5,4.5])


%% Fig. 5b
[~, ord] = sort(median(gid1_unacc_injury_60_post,  2,'omitnan'),'descend');

fig = figure('Color','w');
ax  = axes(fig); hold(ax,'on');

% ---------- medians per region
median_unacc  = [median(gid1_unacc_injury_60_post_risk,  2,'omitnan'),...
                median(gid1_unacc_injury_60_init_risk,  2,'omitnan'),...
                median(gid1_unacc_injury_60_base_risk,  2,'omitnan'),...
                median(gid1_unacc_injury_60_base_init_risk,  2,'omitnan')];
x_lables = {'Impaired (calibrated)','Intact (calibrated)','Impaired (default)','Intact (default)'};
% x_lables = {'Impaired','Intact','Impaired','Intact'};

OUT = plot_unacc_injury_stacked_nature(median_unacc, ...
    city_names, ord, 'ShowTitle', false,...
    'labels',x_lables);

ax.FontName   = 'Helvetica';
ax.FontSize   = 8;
ax.LineWidth  = 0.75;
ax.TickDir    = 'out';
xlim([0.5,4.5])

%% Fig. 5c

addpath("utils");
result_dir = "figures\results";

load(fullfile(result_dir,"all_cell_risk_results_accessibility.mat"),'all_cell_unacc_injury_10_post',...
    'all_cell_unacc_injury_30_post','all_cell_unacc_injury_60_post');

load(fullfile(result_dir,"all_cell_risk_results_accessibility_base.mat"),'all_cell_unacc_injury_10_base',...
    'all_cell_unacc_injury_30_base','all_cell_unacc_injury_60_base', ...
    'cell_informs_uni','all_cell_loss_post');

[gmf_data,gmf_R] = readgeoraster(fullfile(result_dir,"Myanmar_wc_wgs84.tif"));
wealth_values = get_value_from_tif(double(cell_informs_uni(:,[2,3])),gmf_data,gmf_R);


plot_unacc_injury_back2back_grouped_swapped( ...
    wealth_values, ...
    all_cell_unacc_injury_10_base, all_cell_unacc_injury_10_post, ...
    all_cell_unacc_injury_30_base, all_cell_unacc_injury_30_post, ...
    all_cell_unacc_injury_60_base, all_cell_unacc_injury_60_post, ...
    'TierLabels',{'poor','middle','rich'}, 'ThresholdLabels',{'10','30','60'},'UseShare', false);


%% Functions

function OUT = plot_unacc_injury_stacked_nature( ...
    MED, city_names, ord, varargin)
% Nature-style stacked shares of inaccessible injured population.
% - Regions are sorted by their share in IMPAIRED 60 min (median-based).
% - Six stacked columns: Intact 10/30/60, Impaired 10/30/60.
% - Consistent, colorblind-safe palette; minimal, clean styling.
%
% Inputs: each 15x1500 or 1500x15 (rows=regions, cols=draws).
%
% Name-Value:
%   'FontSize'   (default 9)
%   'Colors'     (15x3 RGB; default = Tableau-15 palette below)
%   'ShowTitle'  (true/false; default false)
%   'Title'      (string; used if ShowTitle=true)

ip = inputParser;
ip.addParameter('FontSize',8,@(x)isnumeric(x)&&isscalar(x));
ip.addParameter('Colors',[],@(c)isempty(c)||(isnumeric(c)&&size(c,2)==3));
ip.addParameter('ShowTitle',false,@islogical);
ip.addParameter('Title','Inaccessible injured population — regional shares (median-based)',@(s)ischar(s)||isstring(s));
ip.addParameter('labels','',@(s)iscell(s)||iscellstr(s));


ip.parse(varargin{:});
opt = ip.Results;

% % ---------- ensure orientation (15 rows)
% imp10  = fixrows(imp10); 
% 
% 
% 
% % order bars: Intact 10/30/60, then Impaired 10/30/60
% MED = [m_imp10]; % 15 x 6

% ---------- convert to shares (%) within each bar
PCT = zeros(size(MED));
for j = 1:size(MED,2)
    s = sum(MED(:,j),'omitnan');
    PCT(:,j) = 100 * (MED(:,j) / max(s, eps));
end

% ---------- sort regions by Impaired 60 share (descending)
PCT  = PCT(ord,:);
MED  = MED(ord,:);
if ~isempty(city_names)
    city_names = city_names(:);
    city_names = city_names(ord);
else
    city_names = arrayfun(@(k)sprintf('Region %d',k),1:15,'uni',0).';
end

% ---------- palette (15 distinct, color-blind friendly; Tableau 20 subset)
if isempty(opt.Colors)
    hex = [ ...
        "4E79A7","F28E2B","E15759","76B7B2","59A14F", ...
        "EDC948","B07AA1","FF9DA7","9C755F","BAB0AB", ...
        "1F77B4","2CA02C","D62728","9467BD","17BECF"];
    C = hex2rgb(hex);
else
    C = opt.Colors;
    assert(size(C,1)>=15,'Colors must be 15x3 or more.');
    C = C(1:15,:);
end

% % ---------- plot (Nature-ish styling)


hb = bar(PCT.', 'stacked', 'BarWidth', 0.78);
for k = 1:numel(hb) % one series per region (after transpose)
    hb(k).FaceColor = C(k,:);
    hb(k).EdgeColor = 'none';
end

% axes/labels
ylim([0 100]); yticks(0:20:100);
ylabel('Share of national inaccessible injuries (%)');
xticks(1:numel(opt.labels)); xticklabels(opt.labels);
box 'off'
grid 'off'

% % group divider and caption
% xline(3.5,'--','Color',[0 0 0],'LineWidth',0.6);
% text(2,102,'Intact','HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',opt.FontSize);
% text(5,102,'Impaired','HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',opt.FontSize);

% legend (sorted by Impaired 60 share)
% lg = legend(hb, city_names, 'Location','southoutside', 'NumColumns',2);
% lg.Box = 'off';
% lg.Title.String = 'Regions';

% tidy limits
% xlim([0.4 6.6]);

% ---------- outputs
OUT = struct();
OUT.pct       = PCT;          % 15 x 6 shares (%), sorted by impaired-60
OUT.median    = MED;          % 15 x 6 medians, same order
OUT.order     = ord;          % indices applied to original region order
OUT.colors    = C;            % 15 x 3 RGB used
OUT.bar       = hb;
OUT.labels    = opt.labels;

end

% ---------- helpers
function A = fixrows(A)
    if size(A,1)==15, return; end
    if size(A,2)==15, A = A.'; return; end
    error('Each input must be 15xM or Mx15.');
end

function RGB = hex2rgb(hexStr)
    if isstring(hexStr), hexStr = cellstr(hexStr); end
    n = numel(hexStr); RGB = zeros(n,3);
    for i = 1:n
        h = char(hexStr{i});
        RGB(i,:) = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))] / 255;
    end
end


function OUT = plot_four_pi_bars(samples4xN, yLabelStr)
%PLOT_FOUR_PI_BARS  Bar plot of medians with 95% PI error bars.
%
%   samples4xN : [4 × N] matrix; each row is the results of one series across N simulations
%   yLabelStr  : y-axis label string
%
%   Draw four vertical bars:
%       - bar height = median
%       - error bars = 2.5% and 97.5% percentiles (95% prediction interval, PI)
%   Keep only the y-axis label and ticks; remove everything else
%   (no x-axis ticks, no title, no legend).

    if size(samples4xN,1) ~= 4
        error('Input must be a 4×N matrix.');
    end
    if nargin < 2
        yLabelStr = '';
    end

    % Median & 95% PI
    med = median(samples4xN, 2, 'omitnan');            % 4×1
    prc = prctile(samples4xN.', [2.5 97.5]);           % 2×4 (computed along columns)
    lo  = prc(1,:).';                                  % 4×1
    hi  = prc(2,:).';                                  % 4×1
    errLo = med - lo;
    errHi = hi  - med;

    OUT.med = med;
    OUT.lo = lo;
    OUT.hi = hi;

    % Colours: to avoid clashing with your provided hex colours, use a new set of soft colours
    colors = [ ...
        200 200 200;   % #8DD3C7
        200 200 200;   % #FFFFB3
        200 200 200;   % #BEBADA
        200 200 200];  % #FB8072
    colors = colors / 255;

    % Plot
    figure; hold on;
    x = 1:4;

    b = bar(x, med, 'BarWidth', 0.78);
    b.FaceColor = 'flat';
    for i = 1:4
        b.CData(i,:) = colors(i,:);
    end

    % Error bars (95% PI)
    errorbar(x, med, errLo, errHi, ...
        'k', 'LineStyle', 'none', ...
        'CapSize', 3, 'LineWidth', 1);

    % Keep only the y-axis; remove the rest
    xlim([0.5 4.5]);
    set(gca, ...
        'XTick', [], ...          % do not show x-axis ticks
        'XColor', 'none', ...     % do not draw the x-axis
        'Box', 'off', ...         % remove top and right borders
        'TickDir', 'out', ...
        'LineWidth', 0.75, ...
        'FontName', 'Arial', ...
        'FontSize', 8);

    ylabel(yLabelStr);

    hold off;
end


function fig = plot_unacc_injury_back2back_grouped_swapped( ...
    wealth_values, ...
    u10_base, u10_post, ...
    u30_base, u30_post, ...
    u60_base, u60_post, varargin)

% Back-to-back (top/bottom) grouped bar chart (not stacked) — groups swapped
% x-axis: 10/30/60 minutes (group categories)
% Each group has 3 bars: poor / middle / rich
% Top (positive): calibrated / post
% Bottom (negative): default / base
%
% By default, plots "share": under each threshold + model, poor+middle+rich = 100%

p = inputParser;
addParameter(p,'TierOrder',[1 2 3]);
addParameter(p,'TierLabels',{'poor','middle','rich'});
addParameter(p,'ThresholdLabels',{'10','30','60'});   % x-axis tick labels
addParameter(p,'YLabel','AA injuries with t_h>threshold');
addParameter(p,'UseShare',true); % true = share; false = absolute value
parse(p,varargin{:});

tiers   = p.Results.TierOrder(:)';
tlabels = p.Results.TierLabels;
thrLabs = p.Results.ThresholdLabels;

% --- Vectorise and clean ---
w  = wealth_values(:);
u10b = u10_base(:); u10p = u10_post(:);
u30b = u30_base(:); u30p = u30_post(:);
u60b = u60_base(:); u60p = u60_post(:);

valid = isfinite(w) & isfinite(u10b) & isfinite(u10p) & ...
        isfinite(u30b) & isfinite(u30p) & isfinite(u60b) & isfinite(u60p);

w   = w(valid);
u10b = u10b(valid); u10p = u10p(valid);
u30b = u30b(valid); u30p = u30p(valid);
u60b = u60b(valid); u60p = u60p(valid);

% --- Compute first: rows = poor/middle/rich, columns = 10/30/60 (same as original) ---
S_post = zeros(numel(tiers), 3);
S_base = zeros(numel(tiers), 3);

postSeries = {u10p, u30p, u60p};
baseSeries = {u10b, u30b, u60b};

for j = 1:3
    sp = postSeries{j};
    sb = baseSeries{j};

    if p.Results.UseShare
        denP = max(sum(sp,'omitnan'), eps);
        denB = max(sum(sb,'omitnan'), eps);
    else
        denP = 1; denB = 1;
    end

    for i = 1:numel(tiers)
        ti = tiers(i);
        S_post(i,j) = sum(sp(w==ti),'omitnan') / denP;
        S_base(i,j) = sum(sb(w==ti),'omitnan') / denB;
    end
end

% --- Swap grouping: rows = 10/30/60, columns = poor/middle/rich, so grouped bar uses x=threshold ---
if p.Results.UseShare
    Ypos =  100 * (S_post');    % 3x3 (threshold x poor/middle/rich)
    Yneg = -100 * (S_base');    % 3x3
else
    Ypos =  (S_post');
    Yneg = -(S_base');
end

% --- Plot: x = thresholds, within-group = poor/middle/rich (back-to-back) ---
fig = figure('Color','w');
ax = axes(fig); hold(ax,'on');

x = 1:3;

hPos = bar(ax, x, Ypos, 'grouped', 'BarWidth', 0.78);
hNeg = bar(ax, x, Yneg, 'grouped', 'BarWidth', 0.78);

% Colours correspond to poor/middle/rich (legend too)
cols = [ ...
    0.55 0.10 0.10;   % poor
    0.85 0.35 0.10;   % middle
    0.10 0.25 0.70];  % rich

for i = 1:3
    hPos(i).FaceColor = cols(i,:);
    hPos(i).EdgeColor = 'none';

    hNeg(i).FaceColor = cols(i,:);
    hNeg(i).EdgeColor = 'none';
    hNeg(i).FaceAlpha = 0.35;  % default is shown lighter
end

yline(ax, 0, 'k-', 'LineWidth', 1.2);

set(ax,'XTick',x,'XTickLabel',thrLabs);
xlabel(ax,'threshold (minutes)');
ylabel(ax, p.Results.YLabel);

box(ax,'off');
ax.YGrid = 'on';
ax.XGrid = 'off';

% Show absolute values on y-axis (mirrored top/bottom is clearer)
yt = yticks(ax);
yticklabels(ax, compose('%g', abs(yt)));

legend(ax, tlabels, 'Location','southoutside', 'Orientation','horizontal');

% Annotate top/bottom meaning
xl = xlim(ax); yl = ylim(ax);
text(ax, xl(2), yl(2)*0.85, 'calibrated', 'HorizontalAlignment','right');
text(ax, xl(2), yl(1)*0.85, 'default',    'HorizontalAlignment','right');

end