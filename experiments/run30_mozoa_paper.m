function stats = run30_mozoa_paper(prob, cfg)
%RUN30_MOZOA_PAPER  The exact 5-algorithm comparison used in the MOZOA paper:
%   MOZOA, NSGA-2, MLNSGA-2, SPEA2, and MOEA/D.
%
%   stats = run30_mozoa_paper(prob, cfg)
%
%   This function runs only the five algorithms compared in the paper.
%   MOZOA itself is
%   solveZOA8op.m (the Zebra Optimization Algorithm with the 8-operator
%   neighbourhood, no reinforcement learning), reported here under the
%   name "MOZOA" to match the paper directly.
%
%   NOTE ON BUDGET FAIRNESS: all five algorithms use the same population
%   size and generation/iteration budget by default (cfg.nsga's
%   popSize/numGen), so the comparison is not biased by a larger evaluation
%   budget for any one method.
%
%   INPUTS
%     prob : (optional) instance; defaults to E-n23-k3 ('conflict'), so
%            this can be run with no arguments. For the full 13-instance
%            paper sweep, use run_all_instances_for_paper.m instead, which
%            calls this function once per instance.
%     cfg  : optional struct
%              dataset       'E-n23-k3.vrp' (used only if prob is omitted)
%              numRuns       30 (paper protocol)
%              resume        true
%              referenceAlgo which algorithm's row the C-metric/Contribution/
%                            Wilcoxon columns compare every other algorithm
%                            against (default: 'MOZOA')
%              nsga      solveMOTDGVRP opts, shared base for NSGA-2/MLNSGA-2
%              spea2     solveSPEA2 opts (default: derived from cfg.nsga's
%                        popSize/numGen so the budget matches)
%              moead     solveMOEAD opts (default: derived from cfg.nsga's
%                        popSize/numGen so the budget matches)
%              mozoa      solveZOA8op opts (default: derived from cfg.nsga's
%                        popSize/numGen; opSelectMode/fixedProbs default to
%                        solveZOA8op.m's own defaults, i.e. 'fixed' with the
%                        non-uniform schedule weighted toward O3/O6/O7)
%              verbose   true
%
%   OUTPUT (stats)
%     algos     = {'MOZOA','NSGA-2','MLNSGA-2','SPEA2','MOEA/D'}
%     HV, NDS, CPU : R x 5 matrices
%     bestFront : 1x5 cell, each algorithm's best-run Pareto front (for the
%                 paper's 3-D Pareto-front figures via PLOT_PARETO_FRONTS_3D.m)
%     summary   : per-algorithm best/worst/mean/std/median of HV/NDS/CPU
%     cMetric   : struct with C(MOZOA, X) and C(X, MOZOA) means, X in the
%                 other four algorithms
%     contribution : struct with Contribution(MOZOA vs X) means
%     wilcoxon  : struct with p_MOZOA_vs_<X> for each other algorithm X
%
%   USAGE
%     stats = run30_mozoa_paper;                       % default instance
%     c.numRuns = 30; stats = run30_mozoa_paper([], c);
%
%   For the full paper sweep across all 13 Christofides Set E instances,
%   use run_all_instances_for_paper.m, which is already wired to call this
%   function once per instance and aggregate the results.

if nargin < 2 || isempty(cfg), cfg = struct; end
def = struct('dataset','E-n23-k3.vrp', 'numRuns',30, 'resume',true, ...
    'referenceAlgo', 'MOZOA', ...
    'nsga',   struct('popSize',100,'numGen',100,'verbose',false), ...
    'spea2',  [], 'moead', [], 'mozoa', [], ...
    'verbose', true);
cfg = mergeOptsMozoaPaper(def, cfg);

if isempty(cfg.spea2)
    cfg.spea2 = struct('popSize',cfg.nsga.popSize,'numGen',cfg.nsga.numGen, ...
                        'verbose',false);
end
if isempty(cfg.moead)
    cfg.moead = struct('popSize',cfg.nsga.popSize,'numGen',cfg.nsga.numGen, ...
                        'verbose',false);
end
if isempty(cfg.mozoa)
    cfg.mozoa = struct('popSize',cfg.nsga.popSize,'numIter',cfg.nsga.numGen, ...
                       'verbose',false);
end

if nargin < 1 || isempty(prob)
    if cfg.verbose
        fprintf('No instance given; loading default %s ...\n', cfg.dataset);
    end
    prob = loadCVRPInstance(cfg.dataset, 'conflict', 7);
end

R = cfg.numRuns;
algos = {'MOZOA','NSGA-2','MLNSGA-2','SPEA2','MOEA/D'};
nA = numel(algos);

refIdx = find(strcmp(algos, cfg.referenceAlgo), 1);
if isempty(refIdx)
    error('run30_mozoa_paper:badReference', ...
        'cfg.referenceAlgo = "%s" is not among the algorithms being run (%s).', ...
        cfg.referenceAlgo, strjoin(algos, ', '));
end

fronts = cell(R,nA);
CPU  = zeros(R,nA);

ckpt = sprintf('run30mozoa_checkpoint_%s.mat', regexprep(prob.name,'[^\w]','_'));
startRun = 1;
if cfg.resume && exist(ckpt,'file')
    S = load(ckpt);
    if isfield(S,'fronts') && isfield(S,'completed') && size(S.fronts,1) == R ...
            && size(S.fronts,2) == nA
        fronts = S.fronts; CPU = S.CPU; startRun = S.completed + 1;
        if cfg.verbose
            fprintf('Resuming from checkpoint "%s": %d/%d runs already done.\n', ...
                ckpt, S.completed, R);
        end
    end
end

if cfg.verbose
    fprintf('=== %d runs on %s (n=%d, Q=%d, NV=%d), %d algorithm(s): %s ===\n', ...
        R, prob.name, prob.n, prob.Q, prob.NV, nA, strjoin(algos, ', '));
    fprintf('Reference algorithm for C-metric/Contribution/Wilcoxon: %s\n', cfg.referenceAlgo);
end

for r = startRun:R
    timesThisRun = zeros(1,nA);
    for a = 1:nA
        switch algos{a}
            case 'MOZOA'
                o = cfg.mozoa; o.seed = r;
                res = solveZOA8op(prob, o);
            case 'NSGA-2'
                o = cfg.nsga; o.useML = false; o.seed = r;
                res = solveMOTDGVRP(prob, o);
            case 'MLNSGA-2'
                o = cfg.nsga; o.useML = true; o.seed = r;
                res = solveMOTDGVRP(prob, o);
            case 'SPEA2'
                o = cfg.spea2; o.seed = r;
                res = solveSPEA2(prob, o);
            case 'MOEA/D'
                o = cfg.moead; o.seed = r;
                res = solveMOEAD(prob, o);
        end
        fronts{r,a} = res.paretoObj;
        CPU(r,a) = res.cpuTime;
        timesThisRun(a) = res.cpuTime;
    end

    completed = r; %#ok<NASGU>
    save(ckpt, 'fronts', 'CPU', 'completed', 'R', 'algos', '-v7');

    if cfg.verbose
        parts = arrayfun(@(a) sprintf('%s %.1fs', algos{a}, timesThisRun(a)), ...
            1:nA, 'UniformOutput', false);
        fprintf('  run %2d/%d done (%s)  [saved]\n', r, R, strjoin(parts, ', '));
    end
end

% ---- shared reference point across ALL algorithms/runs ----
pool = [];
for r = 1:R
    for a = 1:nA, pool = [pool; fronts{r,a}]; end %#ok<AGROW>
end
ref = max(pool,[],1)*1.1 + 1;

HV  = zeros(R,nA);
NDS = zeros(R,nA);
for r = 1:R
    for a = 1:nA
        HV(r,a)  = hypervolume(fronts{r,a}, ref);
        NDS(r,a) = size(fronts{r,a},1);
    end
end

stats.algos = algos; stats.ref = ref;
stats.HV = HV; stats.NDS = NDS; stats.CPU = CPU;

% ---- best-run Pareto front per algorithm (for the paper's 3-D PF figures) --
stats.bestFront = cell(1,nA);
stats.bestRunIdx = zeros(1,nA);
for a = 1:nA
    [~, bestR] = max(HV(:,a));
    stats.bestFront{a} = fronts{bestR,a};
    stats.bestRunIdx(a) = bestR;
end

stats.summary = struct();
metrics = {'HV', HV; 'NDS', NDS; 'CPU', CPU};
if cfg.verbose
    fprintf('\n%-10s %-10s %10s %10s %10s %10s %10s\n', ...
        'Metric','Algo','best','worst','mean','std','median');
end
for mIdx = 1:size(metrics,1)
    name = metrics{mIdx,1}; M = metrics{mIdx,2};
    higherBetter = ~strcmp(name,'CPU');
    for a = 1:nA
        col = M(:,a);
        if higherBetter, bestv = max(col); worstv = min(col);
        else,            bestv = min(col); worstv = max(col); end
        s = struct('best',bestv,'worst',worstv,'mean',mean(col), ...
                    'std',std(col),'median',median(col));
        stats.summary.(name)(a) = s;
        if cfg.verbose
            fprintf('%-10s %-10s %10.4g %10.4g %10.4g %10.4g %10.4g\n', ...
                name, algos{a}, bestv, worstv, mean(col), std(col), median(col));
        end
    end
end

% ---- pairwise C-metric & Contribution: MOZOA vs every other algorithm ----
stats.cMetric = struct();
stats.contribution = struct();
refName = algos{refIdx};
refFieldTag = safeFieldNameMozoaPaper(refName);
if cfg.verbose
    fprintf('\n--- Mean coverage C-metric & Contribution over %d runs (reference: %s) ---\n', ...
        R, refName);
end
for a = 1:nA
    if a == refIdx, continue; end
    cZX = zeros(R,1); cXZ = zeros(R,1); contrZX = zeros(R,1);
    for r = 1:R
        Z = fronts{r,refIdx}; X = fronts{r,a};
        cZX(r) = cMetric(Z,X); cXZ(r) = cMetric(X,Z);
        contrZX(r) = contributionMetric(Z,X);
    end
    xTag = safeFieldNameMozoaPaper(algos{a});
    stats.cMetric.(sprintf('C_%s_%s', refFieldTag, xTag)) = mean(cZX);
    stats.cMetric.(sprintf('C_%s_%s', xTag, refFieldTag)) = mean(cXZ);
    stats.contribution.(sprintf('Contrib_%s_vs_%s', refFieldTag, xTag)) = mean(contrZX);
    if cfg.verbose
        fprintf('  C(%s, %-8s) = %.3f   C(%-8s, %s) = %.3f   Contribution = %.3f\n', ...
            refName, algos{a}, mean(cZX), algos{a}, refName, mean(cXZ), mean(contrZX));
    end
end

% ---- Wilcoxon: MOZOA vs every other algorithm on HV ----
stats.wilcoxon = struct();
if exist('signrank','file') == 2
    if cfg.verbose
        fprintf('\n--- Wilcoxon signed-rank on HV (alpha=0.05), %s vs each ---\n', refName);
    end
    for a = 1:nA
        if a == refIdx, continue; end
        xTag = safeFieldNameMozoaPaper(algos{a});
        fname = sprintf('p_%s_vs_%s', refFieldTag, xTag);
        try
            p = signrank(HV(:,refIdx), HV(:,a));
        catch
            p = NaN;
        end
        stats.wilcoxon.(fname) = p;
        if cfg.verbose
            fprintf('  %s vs %-8s: p = %.4g  %s\n', refName, algos{a}, p, sigMozoaPaper(p));
        end
    end
elseif cfg.verbose
    fprintf('\n(Wilcoxon tests skipped: Statistics Toolbox not available.)\n');
end

% ---- HV boxplot-style bar of means with error bars ----
try
    figure('Color','w');
    mu = mean(HV,1); sd = std(HV,0,1);
    bar(mu); hold on;
    errorbar(1:nA, mu, sd, 'k', 'linestyle','none','LineWidth',1.2);
    set(gca,'XTickLabel',algos); ylabel('Hypervolume');
    title(sprintf('%s : mean HV over %d runs (MOZOA paper, 5 algorithms)', prob.name, R)); grid on;
catch
end

finalFile = sprintf('run30mozoa_results_%s.mat', regexprep(prob.name,'[^\w]','_'));
save(finalFile, 'stats', '-v7');
xlsxFile = strrep(finalFile, '.mat', '.xlsx');
try
    exportStatsToExcel(stats, xlsxFile);
catch ME
    if cfg.verbose
        fprintf('(Excel export skipped: %s)\n', ME.message);
    end
end
if exist(ckpt,'file'), delete(ckpt); end

if cfg.verbose
    fprintf('\nResults saved to %s and %s (checkpoint removed).\n', finalFile, xlsxFile);
    fprintf('Done.\n');
end
end

% =========================================================================
function s = sigMozoaPaper(p)
if isnan(p), s = '(n/a)';
elseif p < 0.05, s = '(significant)';
else, s = '(not significant)'; end
end

function s = safeFieldNameMozoaPaper(name)
s = regexprep(name, '[^A-Za-z0-9]', '');
if isempty(s) || ~isletter(s(1)), s = ['A' s]; end
end

function o = mergeOptsMozoaPaper(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
