function allResults = run_all_instances_for_paper(cfg)
%RUN_ALL_INSTANCES_FOR_PAPER  Run run30_mozoa_paper.m on every
%   bundled Christofides Set E instance and aggregate the results into the
%   exact tables/figures the MOZOA paper needs.
%
%   allResults = run_all_instances_for_paper(cfg)
%
%   For each of the 13 bundled E-n*.vrp instances, this:
%     1) builds the instance (with the paper's x100 demand/capacity scaling
%        for E-n51/76/101, matching Nyako et al.'s convention),
%     2) runs run30_mozoa_paper (5 algorithms: MOZOA, NSGA-2, MLNSGA-2,
%        SPEA2, MOEA/D) for cfg.numRuns independent seeds,
%     3) saves that instance's full stats (HV/NDS/CPU matrices, per-algorithm
%        best-run Pareto front for 3-D plotting, C-metric/Contribution/
%        Wilcoxon vs. the reference algorithm) to
%        mozoa_paper_results_<instance>.mat/.xlsx,
%     4) accumulates a per-instance summary row into allResults and saves
%        the running total to mozoa_paper_allresults.mat/.xlsx after every
%        instance, so the whole sweep is safely resumable.
%
%   Per-instance results already on disk are loaded and skipped (not
%   recomputed), so re-running this after an interruption picks up exactly
%   where it left off, and re-running it after adding a new instance file
%   only computes the new one.
%
%   cfg (optional) fields:
%     .instances  cell list of .vrp filenames (default: all 13 bundled
%                 E-n*.vrp, sorted small to large)
%     .numRuns    independent runs per algorithm per instance (default 30)
%     .referenceAlgo  which algorithm is the "proposed method" row for
%                 C-metric/Contribution/Wilcoxon (default 'MOZOA')
%     .scale100   apply the paper's x100 demand/capacity scaling to
%                 E-n51/76/101 (default true)
%     .verbose    default true
%
%   OUTPUT
%     allResults : struct array, one row per instance, with fields
%       name, n, Q, NV, algos, hvMean, hvStd, ndsMean, cpuMean (1 x 5 row
%       vectors, one entry per algorithm in run30_mozoa_paper's algos order:
%       {'MOZOA','NSGA-2','MLNSGA-2','SPEA2','MOEA/D'}), plus the full
%       per-instance stats struct in .stats
%       (so .stats.bestFront{a} gives that algorithm's best-run Pareto
%       front on this instance, for the 3-D PF figures).
%
%   USAGE
%     allResults = run_all_instances_for_paper;             % all 13, 30 runs
%     c.numRuns = 30; allResults = run_all_instances_for_paper(c);
%     c.instances = {'E-n23-k3.vrp','E-n51-k5.vrp'};         % subset only
%     allResults = run_all_instances_for_paper(c);

if nargin < 1 || isempty(cfg), cfg = struct; end

defaultList = { ...
    'E-n13-k4.vrp','E-n22-k4.vrp','E-n23-k3.vrp','E-n30-k3.vrp', ...
    'E-n31-k7.vrp','E-n33-k4.vrp','E-n51-k5.vrp','E-n76-k7.vrp', ...
    'E-n76-k8.vrp','E-n76-k10.vrp','E-n76-k14.vrp','E-n101-k8.vrp', ...
    'E-n101-k14.vrp'};

def = struct('instances',{defaultList}, 'numRuns',30, ...
             'referenceAlgo','MOZOA', 'scale100',true, 'verbose',true);
cfg = mergeOptsPaper(def, cfg);

instances = cfg.instances;
nInst = numel(instances);
allResults = struct('name',{},'n',{},'Q',{},'NV',{},'algos',{}, ...
                    'hvMean',{},'hvStd',{},'ndsMean',{},'cpuMean',{},'stats',{});

if cfg.verbose
    fprintf('##############################################################\n');
    fprintf('# MOZOA paper sweep: %d instance(s) x %d run(s), reference=%s\n', ...
        nInst, cfg.numRuns, cfg.referenceAlgo);
    fprintf('##############################################################\n');
end

sweepStart = tic;

for k = 1:nInst
    fname = instances{k};
    instTag = regexprep(fname, {'\.vrp$','[^\w]'}, {'',''});
    resultFile = sprintf('mozoa_paper_results_%s.mat', instTag);

    if cfg.verbose
        fprintf('\n========== [%d/%d] %s ==========\n', k, nInst, fname);
    end

    if exist(resultFile,'file')
        if cfg.verbose
            fprintf('  found %s -- already done, loading (not re-running the algorithms).\n', resultFile);
        end
        S = load(resultFile);
        stats = S.stats; prob = S.prob;
    else
        prob = loadCVRPInstance(fname, 'conflict', 7);
        if cfg.scale100 && ~isempty(regexp(prob.name,'n(51|76|101)','once'))
            prob = scaleInstance(prob, 100);
            if cfg.verbose, fprintf('  Applied x100 scaling (paper convention).\n'); end
        end

        rc = struct('numRuns',cfg.numRuns, 'referenceAlgo',cfg.referenceAlgo, ...
                    'verbose',cfg.verbose);
        stats = run30_mozoa_paper(prob, rc);

        save(resultFile, 'stats', 'prob', '-v7');
    end

    % Always (re-)export the Excel file, whether stats was just computed or
    % loaded from an existing .mat -- this ensures that if exportStatsToExcel.m
    % is updated (e.g. to add new sheets such as the PF_<algorithm> Pareto-
    % front points), re-running this script on top of already-computed .mat
    % results still regenerates the corresponding .xlsx with the new content,
    % instead of silently keeping a stale Excel file around.
    try
        exportStatsToExcel(stats, strrep(resultFile,'.mat','.xlsx'));
    catch ME
        if cfg.verbose, fprintf('  (Excel export skipped: %s)\n', ME.message); end
    end

    rec.name = prob.name; rec.n = prob.n; rec.Q = prob.Q; rec.NV = prob.NV;
    rec.algos = stats.algos;
    rec.hvMean  = mean(stats.HV,1);  rec.hvStd  = std(stats.HV,0,1);
    rec.ndsMean = mean(stats.NDS,1);
    rec.cpuMean = mean(stats.CPU,1);
    rec.stats = stats;
    allResults(end+1) = rec; %#ok<AGROW>

    save('mozoa_paper_allresults.mat', 'allResults', '-v7');
    try
        exportPaperSummary(allResults, 'mozoa_paper_allresults.xlsx');
    catch ME
        if cfg.verbose, fprintf('  (Summary Excel export skipped: %s)\n', ME.message); end
    end
end

if cfg.verbose
    fprintf('\n##############################################################\n');
    fprintf('# Sweep complete: %d instance(s) in %.1f min\n', ...
        numel(allResults), toc(sweepStart)/60);
    fprintf('# Per-instance files: mozoa_paper_results_<instance>.mat/.xlsx\n');
    fprintf('# Combined summary:   mozoa_paper_allresults.mat/.xlsx\n');
    fprintf('##############################################################\n');
end
end

% =========================================================================
function exportPaperSummary(allResults, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
if isempty(allResults), return; end
algos = allResults(1).algos;
nA = numel(algos);

% ---- HV summary sheet: one row per instance, one column per algorithm ----
header = [{'Instance','n'}, algos];
rows = cell(numel(allResults), numel(header));
for k = 1:numel(allResults)
    r = allResults(k);
    rows(k,:) = [{r.name, r.n}, num2cell(r.hvMean)];
end
writecell([header; rows], xlsxFile, 'Sheet', 'HV_mean');

% ---- NDS summary sheet ----
rows2 = cell(numel(allResults), numel(header));
for k = 1:numel(allResults)
    r = allResults(k);
    rows2(k,:) = [{r.name, r.n}, num2cell(r.ndsMean)];
end
writecell([header; rows2], xlsxFile, 'Sheet', 'NDS_mean');

% ---- CPU summary sheet ----
rows3 = cell(numel(allResults), numel(header));
for k = 1:numel(allResults)
    r = allResults(k);
    rows3(k,:) = [{r.name, r.n}, num2cell(r.cpuMean)];
end
writecell([header; rows3], xlsxFile, 'Sheet', 'CPU_mean');
end

function o = mergeOptsPaper(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
