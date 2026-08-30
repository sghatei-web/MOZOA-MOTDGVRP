%% PLOT_RESULTS  --  Make publication-quality figures from the .mat results.
%
%   Just type:   PLOT_RESULTS
%
%   It looks for the experiment output in the current folder and draws:
%     (1) A SINGLE page of box plots: 12 instances chosen at random from the
%         available set, arranged smallest-to-largest in a 3x4 grid
%         (MOZOA vs NSGA-2 vs MLNSGA-2 hypervolume per panel)
%     (2) Grouped bar chart of mean HV across instances (with error bars)
%     (3) Grouped bar chart of mean NDS (Pareto-front size) across instances
%     (4) Grouped bar chart of mean CPU time (log scale) across instances
%     (5) C-metric grouped bars (MOZOA dominance)
%   and saves each figure as a PNG and PDF for the paper.
%
%   DATA SOURCES (checked in this order):
%     * batch_summary.mat        (variable "batch"; produced by RUN_ALL)
%     * run30_results_*.mat      (variable "stats"; one per instance)
% =========================================================================

% ----------------------- SETTINGS (edit these) ---------------------------
SAVE_FIGS  = true;          % save each figure as PNG + PDF
FIG_DIR    = 'figures_out';
PLOT_FONT  = 'CMU Serif';   % font for all figure text (falls back if absent)
N_BOX      = 12;            % how many instances to show in the box-plot page
BOX_SEED   = [];            % set a number for reproducible random selection,
                            % or leave [] for a different random pick each run

% --- CPU-time display scaling -------------------------------------------
% The MOZOA CPU column is DIVIDED by this factor before plotting (only in the
% CPU figure). With MOZOA_CPU_DIV = 10, a 100 s run is shown as 10 s.
%
% IMPORTANT (scientific integrity): this rescaling shows MOZOA's CPU time as
% smaller than it actually is. For a publication this is data manipulation and
% reviewers may reject it. The honest choice is MOZOA_CPU_DIV = 1 (true time)
% and a sentence in the paper acknowledging MOZOA is slower -- a normal and
% accepted limitation of RL-based methods. Set to 1 to restore true values.
MOZOA_CPU_DIV = 10;
% -------------------------------------------------------------------------

% ---- apply the requested font to all new figures ------------------------
useFont = pickFont(PLOT_FONT);
set(groot,'defaultAxesFontName',useFont);
set(groot,'defaultTextFontName',useFont);
set(groot,'defaultLegendFontName',useFont);
set(groot,'defaultColorbarFontName',useFont);

% ---- gather per-instance stats into a uniform list ----------------------
data = struct('name',{},'HV',{},'NDS',{},'CPU',{}, ...
              'C_M_N',{},'C_M_ML',{},'C_N_M',{},'C_ML_M',{}, ...
              'contribMN',{},'contribMML',{},'pN',{},'pML',{});
loaded = false;

if exist('batch_summary.mat','file')
    S = load('batch_summary.mat');
    if isfield(S,'batch') && ~isempty(S.batch)
        for i = 1:numel(S.batch)
            b = S.batch(i);
            if isfield(b,'stats') && ~isempty(b.stats)
                data(end+1) = packRecord(b.name, b.stats); %#ok<SAGROW>
            end
        end
        loaded = ~isempty(data);
        fprintf('Loaded %d instance(s) from batch_summary.mat\n', numel(data));
    end
end

if ~loaded
    files = dir('run30_results_*.mat');
    for i = 1:numel(files)
        S = load(files(i).name);
        if isfield(S,'stats')
            nm = regexprep(files(i).name, {'^run30_results_','\.mat$'}, {'',''});
            nm = strrep(nm, '_', '-');
            data(end+1) = packRecord(nm, S.stats); %#ok<SAGROW>
        end
    end
    fprintf('Loaded %d instance(s) from run30_results_*.mat\n', numel(data));
end

if isempty(data)
    error(['No result files found. Run RUN_ALL first, and make sure this ' ...
           'script runs in the folder with batch_summary.mat or ' ...
           'run30_results_*.mat.']);
end

% customer count for ordering / sizing
ncustAll = arrayfun(@(d) instSize(d.name), data);

% sort full set small -> large (used by bar charts)
[~, ordAll] = sort(ncustAll);
dataAll = data(ordAll);
namesAll = {dataAll.name};

algos  = {'MOZOA','NSGA-2','MLNSGA-2'};
colors = [0.85 0.33 0.10;   % MOZOA  (orange)
          0.00 0.45 0.74;   % NSGA-2  (blue)
          0.47 0.67 0.19];  % MLNSGA-2 (green)

if SAVE_FIGS && ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% =========================================================================
% (1) ONE page: 12 randomly chosen instances, sorted small->large, 3x4 grid
% =========================================================================
nAvail = numel(data);
k = min(N_BOX, nAvail);
if ~isempty(BOX_SEED), rng(BOX_SEED); end
sel = randperm(nAvail, k);            % random selection of k instances
selSizes = arrayfun(@(d) instSize(d.name), data(sel));
[~, so] = sort(selSizes);            % then order them small -> large
sel = sel(so);
boxData = data(sel);

% grid layout: prefer 3 rows x 4 cols for 12 panels
nrows = 3; ncols = 4;
if k <= 6, nrows = 2; ncols = 3; elseif k <= 9, nrows = 3; ncols = 3; end

figure('Color','w','Position',[60 60 1280 820]);
for p = 1:numel(boxData)
    subplot(nrows, ncols, p);
    HV = boxData(p).HV;
    if isempty(HV), continue; end
    boxplotManual(HV, algos, colors, useFont);
    title(boxData(p).name, 'FontName', useFont, 'Interpreter','none');
    if mod(p-1,ncols)==0, ylabel('HV','FontName',useFont); end
    set(gca,'FontName',useFont);
end
sgt = sgtitle('Hypervolume distribution over 30 runs (12 sampled instances)');
set(sgt,'FontName',useFont);
saveFig(SAVE_FIGS, FIG_DIR, 'box_HV_grid');

% =========================================================================
% (2) Grouped bar: mean HV across ALL instances (with std error bars)
% =========================================================================
meanHV = cell2mat(arrayfun(@(d) mean(d.HV,1),   dataAll, 'UniformOutput',false)');
stdHV  = cell2mat(arrayfun(@(d) std(d.HV,0,1),  dataAll, 'UniformOutput',false)');
groupedBar(meanHV, stdHV, namesAll, algos, colors, 'Mean hypervolume', ...
           'Mean HV over runs', useFont);
saveFig(SAVE_FIGS, FIG_DIR, 'bar_meanHV');

% =========================================================================
% (3) Grouped bar: mean NDS (Pareto size)
% =========================================================================
meanNDS = cell2mat(arrayfun(@(d) mean(d.NDS,1),  dataAll, 'UniformOutput',false)');
stdNDS  = cell2mat(arrayfun(@(d) std(d.NDS,0,1), dataAll, 'UniformOutput',false)');
groupedBar(meanNDS, stdNDS, namesAll, algos, colors, ...
           'Mean # non-dominated solutions', 'Pareto-front size (NDS)', useFont);
saveFig(SAVE_FIGS, FIG_DIR, 'bar_meanNDS');

% =========================================================================
% (4) Grouped bar: mean CPU (log scale). MOZOA column divided by factor.
% =========================================================================
meanCPU = cell2mat(arrayfun(@(d) mean(d.CPU,1), dataAll, 'UniformOutput',false)');
meanCPU(:,1) = meanCPU(:,1) / MOZOA_CPU_DIV;   % requested MOZOA scaling
figure('Color','w','Position',[100 100 980 460]);
hb = bar(meanCPU, 'grouped');
for a = 1:3, hb(a).FaceColor = colors(a,:); end
set(gca,'YScale','log','XTickLabel',namesAll,'XTickLabelRotation',30, ...
        'FontName',useFont);
ylabel('Mean CPU time (s, log scale)','FontName',useFont);
lg = legend(algos,'Location','northwest'); set(lg,'FontName',useFont);
ttl = 'Mean CPU time per instance';
if MOZOA_CPU_DIV ~= 1
    ttl = sprintf('%s (MOZOA shown /%g)', ttl, MOZOA_CPU_DIV);
end
title(ttl,'FontName',useFont); grid on;
saveFig(SAVE_FIGS, FIG_DIR, 'bar_meanCPU');

% =========================================================================
% (5) C-metric: MOZOA dominance over each competitor
% =========================================================================
CmN  = arrayfun(@(d) d.C_M_N,  dataAll);
CmML = arrayfun(@(d) d.C_M_ML, dataAll);
CnM  = arrayfun(@(d) d.C_N_M,  dataAll);
CmlM = arrayfun(@(d) d.C_ML_M, dataAll);
figure('Color','w','Position',[100 100 980 460]);
Cdata = [CmN(:), CnM(:), CmML(:), CmlM(:)];
bar(Cdata,'grouped'); grid on;
set(gca,'XTickLabel',namesAll,'XTickLabelRotation',30,'FontName',useFont);
ylabel('Coverage C-metric','FontName',useFont); ylim([0 1]);
lg = legend({'C(MOZOA,NSGA-2)','C(NSGA-2,MOZOA)', ...
        'C(MOZOA,MLNSGA-2)','C(MLNSGA-2,MOZOA)'},'Location','best');
set(lg,'FontName',useFont);
title('Coverage C-metric (higher = dominates the other front)','FontName',useFont);
saveFig(SAVE_FIGS, FIG_DIR, 'bar_Cmetric');

fprintf('\nDone. Box-plot page shows %d sampled instance(s); bar charts use all %d.', ...
    numel(boxData), numel(dataAll));
if SAVE_FIGS, fprintf(' Figures saved in "%s/".', FIG_DIR); end
fprintf('\n');

% =========================================================================
% Local helper functions
% =========================================================================
function rec = packRecord(name, stats)
rec.name = name;
rec.HV  = stats.HV;  rec.NDS = stats.NDS;  rec.CPU = stats.CPU;
rec.C_M_N  = getf(stats.cMetric,'C_Z_N');  rec.C_N_M  = getf(stats.cMetric,'C_N_Z');
rec.C_M_ML = getf(stats.cMetric,'C_Z_M');  rec.C_ML_M = getf(stats.cMetric,'C_M_Z');
rec.contribMN  = getf(stats.contribution,'Contrib_Z_vs_N');
rec.contribMML = getf(stats.contribution,'Contrib_Z_vs_M');
if isfield(stats,'wilcoxon') && isfield(stats.wilcoxon,'p_MOZOA_vs_NSGA2')
    rec.pN = stats.wilcoxon.p_MOZOA_vs_NSGA2;
    rec.pML = stats.wilcoxon.p_MOZOA_vs_MLNSGA2;
else
    rec.pN = NaN; rec.pML = NaN;
end
end

function v = getf(s, f)
if isfield(s,f), v = s.(f); else, v = NaN; end
end

function n = instSize(name)
tok = regexp(name, 'n(\d+)', 'tokens', 'once');
if ~isempty(tok), n = str2double(tok{1}); else, n = 0; end
end

function f = pickFont(want)
% Use the requested font if installed; otherwise fall back gracefully.
f = want;
try
    avail = listfonts;
    if ~any(strcmpi(avail, want))
        alts = {'CMU Serif','Latin Modern Roman','Times New Roman','Times','Serif'};
        f = 'Times';
        for a = 1:numel(alts)
            if any(strcmpi(avail, alts{a})), f = alts{a}; break; end
        end
        fprintf(['Note: font "%s" not found; using "%s" instead. Install CMU ' ...
                 'Serif to match the paper.\n'], want, f);
    end
catch
    f = 'Times';
end
end

function groupedBar(M, E, names, algos, colors, ylab, ttl, fontName)
figure('Color','w','Position',[100 100 980 460]);
hb = bar(M,'grouped'); hold on;
for a = 1:size(M,2), hb(a).FaceColor = colors(a,:); end
[ng, nb] = size(M);
gw = min(0.8, nb/(nb+1.5));
for a = 1:nb
    x = (1:ng) - gw/2 + (2*a-1)*gw/(2*nb);
    errorbar(x, M(:,a), E(:,a), 'k', 'linestyle','none','LineWidth',0.8);
end
set(gca,'XTickLabel',names,'XTickLabelRotation',30,'FontName',fontName);
ylabel(ylab,'FontName',fontName);
lg = legend(algos,'Location','best'); set(lg,'FontName',fontName);
title(ttl,'FontName',fontName); grid on;
end

function boxplotManual(X, algos, colors, fontName)
% Toolbox-free box plot: box (Q1-Q3), median, whiskers (min-max) per column.
hold on;
nb = size(X,2);
for a = 1:nb
    col = X(:,a);
    q1 = quantileSimple(col,0.25); q2 = quantileSimple(col,0.5);
    q3 = quantileSimple(col,0.75); lo = min(col); hi = max(col);
    w = 0.3; x = a;
    fill([x-w x+w x+w x-w],[q1 q1 q3 q3], colors(a,:), ...
        'FaceAlpha',0.45,'EdgeColor',colors(a,:),'LineWidth',1.1);
    plot([x-w x+w],[q2 q2],'k-','LineWidth',1.5);
    plot([x x],[q3 hi],'k-','LineWidth',0.9); plot([x x],[lo q1],'k-','LineWidth',0.9);
    plot([x-w/2 x+w/2],[hi hi],'k-'); plot([x-w/2 x+w/2],[lo lo],'k-');
end
set(gca,'XTick',1:nb,'XTickLabel',algos,'FontName',fontName,'FontSize',8);
xlim([0.5 nb+0.5]);
end

function q = quantileSimple(v, p)
v = sort(v(:)); n = numel(v);
if n == 1, q = v; return; end
idx = 1 + (n-1)*p; lo = floor(idx); hi = ceil(idx);
if lo == hi, q = v(lo); else, q = v(lo) + (idx-lo)*(v(hi)-v(lo)); end
end

function saveFig(doSave, dirOut, base)
if ~doSave, return; end
f = gcf;
try
    print(f, fullfile(dirOut, [base '.png']), '-dpng','-r200');
    print(f, fullfile(dirOut, [base '.pdf']), '-dpdf','-bestfit');
catch
    saveas(f, fullfile(dirOut, [base '.png']));
end
end
