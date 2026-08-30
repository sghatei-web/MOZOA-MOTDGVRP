%% PLOT_PAPER_FIGURES_FROM_EXCEL.m
%  Reads the Excel files produced by run_all_instances_for_paper.m /
%  run30_mozoa_paper.m and draws the figures needed for the MOZOA paper:
%    1) 3-D Pareto-front scatter plots (distance/time/fuel), one per
%       instance, all algorithms overlaid
%    2) Box plots of hypervolume (HV) per algorithm, one panel per instance,
%       arranged in a grid across all instances
%    3) Grouped bar charts of mean HV, mean NDS, and mean CPU across all
%       instances
%
%  INPUTS EXPECTED (produced by MATLAB, no manual editing needed):
%    mozoa_paper_allresults.xlsx        -- sheets HV_mean, NDS_mean, CPU_mean
%                                          (one row per instance, one column
%                                          per algorithm)
%    mozoa_paper_results_<instance>.xlsx -- per instance, sheets Raw_HV,
%                                          Raw_NDS, Raw_CPU (one row per run)
%                                          and PF_<Algorithm> (one row per
%                                          Pareto point: Distance/Time/Fuel)
%
%  USAGE: place this script in the same folder as the Excel files and run
%    PLOT_PAPER_FIGURES_FROM_EXCEL
%  All figures are saved as PNG and PDF under figures_paper/.

% ----------------------- SETTINGS (edit these) ---------------------------
SAVE_FIGS  = true;
FIG_DIR    = 'figures_paper';
PLOT_FONT  = 'CMU Serif';      % falls back to Times if not installed
N_BOX      = 12;               % how many instances to show in the box-plot grid
BOX_SEED   = [];                % [] = different random 12-of-13 pick each run
PF_INSTANCES_TO_PLOT = 6;       % how many instances get a 3-D PF figure
% -------------------------------------------------------------------------

useFont = pickFontPaper(PLOT_FONT);
set(groot,'defaultAxesFontName',useFont);
set(groot,'defaultTextFontName',useFont);
set(groot,'defaultLegendFontName',useFont);

if SAVE_FIGS && ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% ---- discover per-instance result files ----------------------------------
files = dir('mozoa_paper_results_*.xlsx');
if isempty(files)
    error(['No mozoa_paper_results_*.xlsx files found in the current folder. ' ...
           'Run run_all_instances_for_paper first.']);
end

instTags  = cell(numel(files),1);
instSizes = zeros(numel(files),1);
for i = 1:numel(files)
    instTags{i} = regexprep(files(i).name, {'^mozoa_paper_results_','\.xlsx$'}, {'',''});
    tok = regexp(instTags{i}, 'n(\d+)', 'tokens', 'once');
    if ~isempty(tok), instSizes(i) = str2double(tok{1}); end
end
[~, ord] = sort(instSizes);
files = files(ord); instTags = instTags(ord); instSizes = instSizes(ord); %#ok<NASGU>

algos  = {'MOZOA','NSGA-2','MLNSGA-2','SPEA2','MOEA/D'};
nA = numel(algos);
colors = [0.85 0.33 0.10;   % MOZOA   -- orange (proposed method)
          0.00 0.45 0.74;   % NSGA-2  -- blue
          0.47 0.67 0.19;   % MLNSGA-2-- green
          0.49 0.18 0.56;   % SPEA2   -- purple
          0.93 0.69 0.13];  % MOEA/D  -- gold
markers = {'o','s','^','d','v'};

% =========================================================================
% (1) 3-D Pareto-front figures
% =========================================================================
nInst = numel(files);
kPF = min(PF_INSTANCES_TO_PLOT, nInst);
selIdxPF = unique(round(linspace(1, nInst, kPF)));

for ii = 1:numel(selIdxPF)
    idx = selIdxPF(ii);
    xlsxFile = fullfile(files(idx).folder, files(idx).name);
    figure('Color','w','Position',[80 80 780 640]); hold on;
    legendEntries = {};
    for a = 1:nA
        sheetName = ['PF_' safeSheetTagPaper(algos{a})];
        try
            T = readcell(xlsxFile, 'Sheet', sheetName);
        catch
            continue;  % this algorithm's sheet may be absent/empty
        end
        if size(T,1) < 2, continue; end
        F = cell2mat(T(2:end,1:3));
        scatter3(F(:,1), F(:,2), F(:,3), 46, colors(a,:), markers{a}, ...
            'filled', 'MarkerEdgeColor','k', 'LineWidth',0.4, 'MarkerFaceAlpha',0.85);
        legendEntries{end+1} = algos{a}; %#ok<SAGROW>
    end
    xlabel('Distance','FontName',useFont); ylabel('Time','FontName',useFont);
    zlabel('Fuel','FontName',useFont);
    title(sprintf('Pareto fronts on %s', instTags{idx}), 'FontName', useFont, 'Interpreter','none');
    legend(legendEntries, 'Location','eastoutside', 'FontName', useFont);
    grid on; view(-37.5, 24); box on;
    if SAVE_FIGS, saveFigPaper(FIG_DIR, sprintf('PF3D_%s', instTags{idx})); end
end

% =========================================================================
% (2) HV box plots: 12-of-13 instances, random selection, sorted small->large
% =========================================================================
kBox = min(N_BOX, nInst);
if ~isempty(BOX_SEED), rng(BOX_SEED); end
selBox = randperm(nInst, kBox);
[~, so] = sort(instSizes(selBox));
selBox = selBox(so);

nrows = 3; ncols = 4;
if kBox <= 6, nrows = 2; ncols = 3; elseif kBox <= 9, nrows = 3; ncols = 3; end

figure('Color','w','Position',[60 60 1280 820]);
for p = 1:numel(selBox)
    idx = selBox(p);
    xlsxFile = fullfile(files(idx).folder, files(idx).name);
    T = readcell(xlsxFile, 'Sheet', 'Raw_HV');
    hdr = T(1,2:end); HV = cell2mat(T(2:end,2:end));
    % align columns to the canonical algos order (defensive against reordering)
    colIdx = zeros(1,nA);
    for a = 1:nA
        f = find(strcmp(hdr, algos{a}), 1);
        if ~isempty(f), colIdx(a) = f; end
    end
    subplot(nrows, ncols, p);
    boxplotManualPaper(HV, colIdx, algos, colors, useFont);
    title(instTags{idx}, 'FontName', useFont, 'Interpreter','none', 'FontSize', 9);
    if mod(p-1,ncols)==0, ylabel('HV','FontName',useFont); end
end
sgt = sgtitle('Hypervolume distribution across instances');
set(sgt,'FontName',useFont);
if SAVE_FIGS, saveFigPaper(FIG_DIR, 'box_HV_grid'); end

% =========================================================================
% (3) Grouped bar charts: mean HV, NDS, CPU across ALL instances
% =========================================================================
summaryFile = 'mozoa_paper_allresults.xlsx';
if exist(summaryFile,'file')
    plotBarSheet(summaryFile, 'HV_mean', 'Mean hypervolume', algos, colors, ...
        useFont, FIG_DIR, 'bar_meanHV', SAVE_FIGS, false);
    plotBarSheet(summaryFile, 'NDS_mean', 'Mean # non-dominated solutions', ...
        algos, colors, useFont, FIG_DIR, 'bar_meanNDS', SAVE_FIGS, false);
    plotBarSheet(summaryFile, 'CPU_mean', 'Mean CPU time (s)', algos, colors, ...
        useFont, FIG_DIR, 'bar_meanCPU', SAVE_FIGS, true);
else
    warning('%s not found; skipping the cross-instance bar charts.', summaryFile);
end

fprintf('\nDone. Figures saved in "%s/".\n', FIG_DIR);

% =========================================================================
% Local helper functions
% =========================================================================
function plotBarSheet(xlsxFile, sheetName, ylab, algos, colors, useFont, ...
                       figDir, baseName, saveFigs, logScale)
T = readcell(xlsxFile, 'Sheet', sheetName);
hdr = T(1,:); body = T(2:end,:);
instNames = body(:,1);
M = cell2mat(body(:,3:end));  % columns: Instance, n, algo1..algoN
colHdr = hdr(3:end);
colIdx = zeros(1,numel(algos));
for a = 1:numel(algos)
    f = find(strcmp(colHdr, algos{a}), 1);
    if ~isempty(f), colIdx(a) = f; end
end
M = M(:, colIdx(colIdx>0));
usedAlgos = algos(colIdx>0);
usedColors = colors(colIdx>0,:);

figure('Color','w','Position',[100 100 1000 460]);
hb = bar(M, 'grouped');
for a = 1:size(M,2), hb(a).FaceColor = usedColors(a,:); end
set(gca,'XTickLabel',instNames,'XTickLabelRotation',30,'FontName',useFont);
if logScale, set(gca,'YScale','log'); end
ylabel(ylab,'FontName',useFont);
lg = legend(usedAlgos,'Location','best'); set(lg,'FontName',useFont);
title(ylab,'FontName',useFont); grid on;
if saveFigs, saveFigPaper(figDir, baseName); end
end

function boxplotManualPaper(HV, colIdx, algos, colors, fontName)
hold on;
nb = numel(algos);
for a = 1:nb
    if colIdx(a) == 0, continue; end
    col = HV(:, colIdx(a));
    q1 = quantileSimplePaper(col,0.25); q2 = quantileSimplePaper(col,0.5);
    q3 = quantileSimplePaper(col,0.75); lo = min(col); hi = max(col);
    w = 0.3; x = a;
    fill([x-w x+w x+w x-w],[q1 q1 q3 q3], colors(a,:), ...
        'FaceAlpha',0.45,'EdgeColor',colors(a,:),'LineWidth',1.1);
    plot([x-w x+w],[q2 q2],'k-','LineWidth',1.5);
    plot([x x],[q3 hi],'k-','LineWidth',0.9); plot([x x],[lo q1],'k-','LineWidth',0.9);
    plot([x-w/2 x+w/2],[hi hi],'k-'); plot([x-w/2 x+w/2],[lo lo],'k-');
end
set(gca,'XTick',1:nb,'XTickLabel',algos,'FontName',fontName,'FontSize',7, ...
    'XTickLabelRotation',30);
xlim([0.5 nb+0.5]);
end

function q = quantileSimplePaper(v, p)
v = sort(v(:)); n = numel(v);
if n == 1, q = v; return; end
idx = 1 + (n-1)*p; lo = floor(idx); hi = ceil(idx);
if lo == hi, q = v(lo); else, q = v(lo) + (idx-lo)*(v(hi)-v(lo)); end
end

function f = pickFontPaper(want)
f = want;
try
    avail = listfonts;
    if ~any(strcmpi(avail, want))
        alts = {'CMU Serif','Latin Modern Roman','Times New Roman','Times'};
        f = 'Times';
        for a = 1:numel(alts)
            if any(strcmpi(avail, alts{a})), f = alts{a}; break; end
        end
    end
catch
    f = 'Times';
end
end

function s = safeSheetTagPaper(name)
s = regexprep(name, '[\\/\?\*\[\]:]', '');
if numel(s) > 27, s = s(1:27); end
end

function saveFigPaper(dirOut, base)
f = gcf;
try
    print(f, fullfile(dirOut, [base '.png']), '-dpng', '-r220');
    print(f, fullfile(dirOut, [base '.pdf']), '-dpdf', '-bestfit');
catch
    saveas(f, fullfile(dirOut, [base '.png']));
end
end
