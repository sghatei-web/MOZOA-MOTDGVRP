%% PLOT_PARETO_FRONTS_3D.m
%  Draws 3-D Pareto-front comparison figures (distance / time / fuel) for the
%  MOZOA paper, one figure per selected instance, with every algorithm's
%  best-run front plotted as a distinct marker/color on the same axes --
%  the 3-objective analogue of Fig. 6/7 in Li et al. (2023) "Deep
%  reinforcement learning for multi-objective combinatorial optimization"
%  (which plots 2-D PFs since that paper studies bi-objective TSP; this
%  problem is tri-objective, so distance/time/fuel are all shown on one 3-D
%  axes per instance instead of two side-by-side 2-D panels).
%
%  PREREQUISITE: run run_all_instances_for_paper first, so that
%  mozoa_paper_results_<instance>.mat files exist in the current folder
%  (each contains stats.bestFront{a}, the best-run Pareto front per
%  algorithm on that instance).
%
%  USAGE:  just type   PLOT_PARETO_FRONTS_3D
%  Customise which instances to plot and how many (paper asked for 4 or 6)
%  via the settings below.

% ----------------------- SETTINGS (edit these) ---------------------------
N_INSTANCES_TO_PLOT = 6;      % paper asked for 4 or 6 instances
INSTANCE_SELECTION  = 'spread'; % 'spread' picks N_INSTANCES_TO_PLOT instances
                                 % evenly spaced from small to large; or give
                                 % an explicit cell list of instance tags
                                 % instead, e.g. {'E_n23_k3','E_n51_k5', ...}
SAVE_FIGS = true;
FIG_DIR   = 'figures_pf3d';
PLOT_FONT = 'CMU Serif';       % falls back to Times if not installed
% -------------------------------------------------------------------------

useFont = pickFontPF(PLOT_FONT);
set(groot,'defaultAxesFontName',useFont);
set(groot,'defaultTextFontName',useFont);
set(groot,'defaultLegendFontName',useFont);

files = dir('mozoa_paper_results_*.mat');
if isempty(files)
    error(['No mozoa_paper_results_*.mat files found. Run ' ...
           'run_all_instances_for_paper first.']);
end

% ---- load every instance's stats, keep instance name + customer count ---
data = struct('tag',{},'name',{},'n',{},'stats',{});
for i = 1:numel(files)
    S = load(files(i).name);
    tag = regexprep(files(i).name, {'^mozoa_paper_results_','\.mat$'}, {'',''});
    data(end+1) = struct('tag',tag,'name',S.prob.name,'n',S.prob.n,'stats',S.stats); %#ok<SAGROW>
end

[~, ord] = sort([data.n]);
data = data(ord);

if ischar(INSTANCE_SELECTION) && strcmpi(INSTANCE_SELECTION,'spread')
    nAvail = numel(data);
    k = min(N_INSTANCES_TO_PLOT, nAvail);
    if k >= nAvail
        selIdx = 1:nAvail;
    else
        selIdx = round(linspace(1, nAvail, k));
        selIdx = unique(selIdx);
    end
else
    selIdx = find(ismember({data.tag}, INSTANCE_SELECTION));
end
selData = data(selIdx);

algos = selData(1).stats.algos;
nA = numel(algos);
colors = [0.85 0.33 0.10;   % ZOA-8op / proposed -- orange (set below)
          0.00 0.45 0.74;   % blue
          0.47 0.67 0.19;   % green
          0.49 0.18 0.56;   % purple
          0.93 0.69 0.13;   % yellow/gold
          0.30 0.30 0.30;   % dark grey
          0.64 0.08 0.18];  % dark red
markers = {'o','s','^','d','v','p','h'};

% put the reference/proposed algorithm ('ZOA-8op' by default) first in the
% legend and give it the most visually prominent colour+marker
propIdx = find(strcmp(algos,'ZOA-8op'), 1);
if isempty(propIdx), propIdx = 1; end
order = [propIdx, setdiff(1:nA, propIdx)];

if SAVE_FIGS && ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

for i = 1:numel(selData)
    d = selData(i);
    figure('Color','w','Position',[80 80 780 640]);
    hold on;
    legendEntries = cell(1,nA);
    for jj = 1:nA
        a = order(jj);
        F = d.stats.bestFront{a};
        if isempty(F), continue; end
        scatter3(F(:,1), F(:,2), F(:,3), 46, colors(jj,:), markers{jj}, ...
            'filled', 'MarkerEdgeColor','k', 'LineWidth',0.4, 'MarkerFaceAlpha',0.85);
        legendEntries{jj} = algos{a};
    end
    xlabel('Distance', 'FontName', useFont);
    ylabel('Time', 'FontName', useFont);
    zlabel('Fuel', 'FontName', useFont);
    title(sprintf('Pareto fronts on %s (n=%d)', d.name, d.n), 'FontName', useFont);
    legend(legendEntries(~cellfun(@isempty,legendEntries)), 'Location','eastoutside', ...
        'FontName', useFont);
    grid on; view(-37.5, 24); box on;
    hold off;

    if SAVE_FIGS
        base = sprintf('PF3D_%s', d.tag);
        try
            print(gcf, fullfile(FIG_DIR, [base '.png']), '-dpng', '-r220');
            print(gcf, fullfile(FIG_DIR, [base '.pdf']), '-dpdf', '-bestfit');
        catch
            saveas(gcf, fullfile(FIG_DIR, [base '.png']));
        end
    end
end

fprintf('Plotted %d instance(s): %s\n', numel(selData), strjoin({selData.tag}, ', '));
if SAVE_FIGS, fprintf('Figures saved in "%s/".\n', FIG_DIR); end

% =========================================================================
function f = pickFontPF(want)
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
