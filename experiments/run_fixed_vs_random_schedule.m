function results = run_fixed_vs_random_schedule(cfg)
%RUN_FIXED_VS_RANDOM_SCHEDULE  Head-to-head comparison of MOZOA's fixed
%   operator-selection schedule (Eq. 40 of the paper) against uniform
%   random operator selection over the same eight-operator set, across
%   all thirteen Christofides Set E instances.
%
%   results = run_fixed_vs_random_schedule(cfg)
%
%   For each instance, runs MOZOA with opSelectMode='fixed' and again with
%   opSelectMode='random', for the same number of independent runs and the
%   same population/iteration budget, then reports mean/std hypervolume
%   per schedule, a paired Wilcoxon signed-rank test (fixed vs random,
%   matched by seed), and which schedule wins on mean HV.
%
%   INPUTS (all optional; cfg is a struct with any of these fields)
%     cfg.instances   cell array of .vrp filenames (default: all 13 bundled)
%     cfg.numRuns     independent runs per schedule per instance (default 30)
%     cfg.popSize     population size (default 100, matching the paper)
%     cfg.numIter     iteration budget (default 100, matching the paper)
%     cfg.verbose     default true
%
%   OUTPUT
%     results : struct array, one row per instance, with fields
%       name, n, hvMeanFixed, hvStdFixed, hvMeanRandom, hvStdRandom,
%       cpuMeanFixed, cpuMeanRandom, wilcoxonP, winner, significant
%     Also saved to fixed_vs_random_results.mat/.xlsx.
%
%   USAGE
%     run_fixed_vs_random_schedule;                    % full: all 13 instances
%     c.numRuns = 30; c.popSize=100; c.numIter=100;
%     run_fixed_vs_random_schedule(c);                  % quick smoke test

if nargin < 1 || isempty(cfg), cfg = struct; end
defaultList = { ...
    'E-n13-k4.vrp','E-n22-k4.vrp','E-n23-k3.vrp','E-n30-k3.vrp', ...
    'E-n31-k7.vrp','E-n33-k4.vrp','E-n51-k5.vrp','E-n76-k7.vrp', ...
    'E-n76-k8.vrp','E-n76-k10.vrp','E-n76-k14.vrp','E-n101-k8.vrp', ...
    'E-n101-k14.vrp'};
def = struct('instances',{defaultList}, 'numRuns',30, 'popSize',100, ...
             'numIter',100, 'verbose',true);
cfg = mergeOptsFixRand(def, cfg);

results = struct('name',{},'n',{},'hvMeanFixed',{},'hvStdFixed',{}, ...
                  'hvMeanRandom',{},'hvStdRandom',{},'cpuMeanFixed',{}, ...
                  'cpuMeanRandom',{},'wilcoxonP',{},'winner',{},'significant',{});

if cfg.verbose
    fprintf('##############################################################\n');
    fprintf('# Fixed vs. random operator schedule: %d instance(s), N=%d runs each\n', ...
        numel(cfg.instances), cfg.numRuns);
    fprintf('##############################################################\n');
end

sweepStart = tic;

for k = 1:numel(cfg.instances)
    fname = cfg.instances{k};
    prob = loadCVRPInstance(fname, 'conflict', 7);
    if ~isempty(regexp(prob.name, 'n(51|76|101)', 'once'))
        prob = scaleInstance(prob, 100);
    end

    if cfg.verbose
        fprintf('\n========== [%d/%d] %s (n=%d) ==========\n', k, numel(cfg.instances), prob.name, prob.n);
    end

    [hvFixed, cpuFixed, frontsFixed] = runScheduleFixRand(prob, 'fixed', cfg.numRuns, cfg.popSize, cfg.numIter);
    [hvRandomRaw, cpuRandom, frontsRandom] = runScheduleFixRand(prob, 'random', cfg.numRuns, cfg.popSize, cfg.numIter); %#ok<ASGLU>

    % shared reference point across BOTH schedules' fronts for a fair
    % paired hypervolume comparison
    pool = [];
    for r = 1:cfg.numRuns
        pool = [pool; frontsFixed{r}; frontsRandom{r}]; %#ok<AGROW>
    end
    ref = max(pool,[],1)*1.1 + 1;

    hvF = zeros(cfg.numRuns,1); hvR = zeros(cfg.numRuns,1);
    for r = 1:cfg.numRuns
        hvF(r) = hypervolume(frontsFixed{r}, ref);
        hvR(r) = hypervolume(frontsRandom{r}, ref);
    end

    if exist('signrank','file') == 2
        try
            p = signrank(hvF, hvR);
        catch
            p = NaN;
        end
    else
        p = NaN;
    end

    if mean(hvF) > mean(hvR), winner = 'fixed'; else, winner = 'random'; end
    sig = ~isnan(p) && p < 0.05;

    if cfg.verbose
        fprintf('  fixed : HV mean=%.4g (std=%.4g), CPU mean=%.2fs\n', mean(hvF), std(hvF), mean(cpuFixed));
        fprintf('  random: HV mean=%.4g (std=%.4g), CPU mean=%.2fs\n', mean(hvR), std(hvR), mean(cpuRandom));
        fprintf('  Wilcoxon (fixed vs random): p=%.4g (%s, favours %s)\n', ...
            p, tern(sig,'significant','not significant'), winner);
    end

    rec = struct('name',prob.name, 'n',prob.n, ...
        'hvMeanFixed',mean(hvF), 'hvStdFixed',std(hvF), ...
        'hvMeanRandom',mean(hvR), 'hvStdRandom',std(hvR), ...
        'cpuMeanFixed',mean(cpuFixed), 'cpuMeanRandom',mean(cpuRandom), ...
        'wilcoxonP',p, 'winner',winner, 'significant',sig);
    results(end+1) = rec; %#ok<AGROW>

    save('fixed_vs_random_results.mat', 'results', '-v7');
    try
        exportFixRandToExcel(results, 'fixed_vs_random_results.xlsx');
    catch ME
        if cfg.verbose, fprintf('(Excel export skipped: %s)\n', ME.message); end
    end
end

nFixedWins = sum(strcmp({results.winner}, 'fixed'));
nSig = sum([results.significant]);
if cfg.verbose
    fprintf('\n##############################################################\n');
    fprintf('# Sweep complete in %.1f min\n', toc(sweepStart)/60);
    fprintf('# Fixed schedule wins on mean HV: %d/%d instances\n', nFixedWins, numel(results));
    fprintf('# Statistically significant differences (uncorrected alpha=0.05): %d/%d instances\n', ...
        nSig, numel(results));
    fprintf('##############################################################\n');
end
end

% =========================================================================
function [hv, cpu, fronts] = runScheduleFixRand(prob, mode, numRuns, popSize, numIter)
hv = zeros(numRuns,1);
cpu = zeros(numRuns,1);
fronts = cell(numRuns,1);
for r = 1:numRuns
    opts = struct('popSize',popSize,'numIter',numIter,'opSelectMode',mode, ...
                  'verbose',false,'seed',r);
    res = solveZOA8op(prob, opts);
    fronts{r} = res.paretoObj;
    cpu(r) = res.cpuTime;
end
end

function s = tern(cond, a, b)
if cond, s = a; else, s = b; end
end

function exportFixRandToExcel(results, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
header = {'Instance','n','HV_mean_fixed','HV_std_fixed','HV_mean_random', ...
          'HV_std_random','CPU_mean_fixed','CPU_mean_random','Wilcoxon_p', ...
          'Winner_by_mean_HV','Significant_at_0.05'};
rows = cell(numel(results), numel(header));
for i = 1:numel(results)
    r = results(i);
    rows(i,:) = {r.name, r.n, r.hvMeanFixed, r.hvStdFixed, r.hvMeanRandom, ...
                 r.hvStdRandom, r.cpuMeanFixed, r.cpuMeanRandom, r.wilcoxonP, ...
                 r.winner, r.significant};
end
writecell([header; rows], xlsxFile, 'Sheet', 'Summary');
end

function o = mergeOptsFixRand(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
