function allResults = run_large_scale_evaluation(cfg)
%RUN_LARGE_SCALE_EVALUATION  Scalability evaluation on larger benchmark
%   instances (Reviewer request: demonstrate how solution quality and
%   runtime scale beyond the 100-customer ceiling of the Christofides
%   Set E benchmark).
%
%   ===================================================================
%   STEP 0 (do this once): DOWNLOAD THE LARGE INSTANCES
%   ===================================================================
%   This function does NOT download instances itself. Fetch a set of
%   larger CVRP instances first -- the Golden et al. benchmark (20
%   instances, 200-480 customers, the standard "large-scale VRP"
%   benchmark used across the VRP literature) is the natural choice.
%   Download the files from
%   https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-instances/
%   ("Golden et al." benchmark section) and place them in a folder, e.g.
%   "instances_large/", next to this script.
%
%   IMPORTANT FORMAT NOTE: the Golden benchmark files from that site
%   (kelly01.txt..kelly20.txt, or all 20 instances back-to-back in
%   kelly-pbs.txt) use a DIFFERENT, non-CVRPLIB-standard text format
%   (no "NAME:"/"DIMENSION:"/"NODE_COORD_SECTION" keywords -- just
%   "[k] with N Customers", "Vehicle capacity = ...", then
%   "Point x y q" rows). readCVRP.m/loadCVRPInstance.m do NOT understand
%   this format (renaming the file to .vrp does not help -- the content
%   itself is different), so this script uses readGolden.m/
%   loadGoldenInstance.m instead, which parse this format directly and
%   still produce a fully compatible MOTDGVRP instance via buildInstance.m.
%   Keep the files as .txt; do not rename them to .vrp.
%
%   ===================================================================
%   STEP 1: run this function
%   ===================================================================
%     run_large_scale_evaluation;                          % default folder
%     c.instancesDir = 'instances_large'; c.numRuns = 30;
%     run_large_scale_evaluation(c);
%     c.numRuns = 30; c.popSize = 100; c.numIter = 100;
%     run_large_scale_evaluation(c);                        % quick smoke test
%
%   INPUTS (all optional; cfg is a struct with any of these fields)
%     cfg.instancesDir  folder containing the downloaded .txt files
%                       (default: 'instances_large')
%     cfg.numRuns       independent runs per algorithm per instance
%                       (default 30, matching the paper protocol)
%     cfg.popSize       population size (default 100, matching the paper)
%     cfg.numIter       iteration/generation budget (default 100)
%     cfg.referenceAlgo which algorithm is the "proposed method" row for
%                       C-metric/Contribution/Wilcoxon (default 'MOZOA')
%     cfg.verbose       default true
%
%   OUTPUT
%     allResults : struct array, one row per instance, with the same
%       fields as run_all_instances_for_paper.m's output (name, n, Q, NV,
%       algos, hvMean, hvStd, ndsMean, cpuMean, stats).
%     Also saves large_scale_results_<instance>.mat/.xlsx (per instance)
%     and large_scale_allresults.mat/.xlsx (cross-instance summary).

if nargin < 1 || isempty(cfg), cfg = struct; end
def = struct('instancesDir','instances_large', 'numRuns',30, ...
    'popSize',100, 'numIter',100, 'referenceAlgo','MOZOA', 'verbose',true);
cfg = mergeOptsLarge(def, cfg);

files = dir(fullfile(cfg.instancesDir, '*.txt'));
if isempty(files)
    error(['No .txt files found in "%s". See the header comment of this ' ...
           'function for how to download large Golden-benchmark instances first.'], ...
           cfg.instancesDir);
end

% ---- warn if BOTH single-instance files (kelly01.txt..kelly20.txt) AND
%      a combined multi-instance file (kelly-pbs.txt) are present, since
%      that would process the same 20 instances twice (once from each
%      source) -- wasteful and confusing when comparing results ----
nSingleFiles = sum(arrayfun(@(f) detectGoldenBlockCount(fullfile(f.folder,f.name))==1, files));
nMultiFiles  = numel(files) - nSingleFiles;
if nSingleFiles > 0 && nMultiFiles > 0
    warning('run_large_scale_evaluation:possibleDuplicates', ...
        ['Found %d single-instance file(s) AND %d multi-instance file(s) in "%s". ' ...
         'If these overlap (e.g. kelly01.txt and kelly-pbs.txt''s block [1] are the ' ...
         'same instance), you will process the same instances twice under different ' ...
         'tags. Consider keeping only one source (either the individual kellyNN.txt ' ...
         'files, or the combined kelly-pbs.txt) to avoid duplicate work.'], ...
        nSingleFiles, nMultiFiles, cfg.instancesDir);
end

% ---- expand every file into one work item per instance block ----------
% Some downloaded files (e.g. kelly01.txt..kelly20.txt) contain exactly one
% instance each; others (e.g. kelly-pbs.txt) contain all 20 back-to-back.
% detectGoldenBlockCount inspects each file once (cheap: just counts header
% lines) so both kinds of files work automatically with no extra options.
workItems = struct('fname',{},'fpath',{},'blockIdx',{},'tag',{});
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    nBlocks = detectGoldenBlockCount(fpath);
    baseTag = regexprep(files(k).name, {'\.txt$','[^\w]'}, {'',''});
    if nBlocks <= 1
        workItems(end+1) = struct('fname',files(k).name, 'fpath',fpath, ...
            'blockIdx',[], 'tag',baseTag); %#ok<AGROW>
    else
        for b = 1:nBlocks
            workItems(end+1) = struct('fname',files(k).name, 'fpath',fpath, ...
                'blockIdx',b, 'tag',sprintf('%s_block%02d', baseTag, b)); %#ok<AGROW>
        end
        if cfg.verbose
            fprintf('  %s contains %d instances -- will process each as a separate item.\n', ...
                files(k).name, nBlocks);
        end
    end
end

allResults = struct('name',{},'n',{},'Q',{},'NV',{},'algos',{}, ...
                    'hvMean',{},'hvStd',{},'ndsMean',{},'cpuMean',{},'stats',{});

if cfg.verbose
    fprintf('##############################################################\n');
    fprintf('# Large-scale scalability sweep: %d instance(s), N=%d runs, pop=%d, iter=%d\n', ...
        numel(workItems), cfg.numRuns, cfg.popSize, cfg.numIter);
    fprintf('##############################################################\n');
end

sweepStart = tic;

for k = 1:numel(workItems)
    w = workItems(k);
    resultFile = sprintf('large_scale_results_%s.mat', w.tag);

    if cfg.verbose
        if isempty(w.blockIdx)
            fprintf('\n========== [%d/%d] %s ==========\n', k, numel(workItems), w.fname);
        else
            fprintf('\n========== [%d/%d] %s (block %d) ==========\n', k, numel(workItems), w.fname, w.blockIdx);
        end
    end

    if exist(resultFile,'file')
        if cfg.verbose
            fprintf('  found %s -- already done, loading (not re-running).\n', resultFile);
        end
        S = load(resultFile);
        stats = S.stats; prob = S.prob;
    else
        tLoad = tic;
        prob = loadGoldenInstance(w.fpath, w.blockIdx, 'conflict', 7);
        if cfg.verbose
            fprintf('  Loaded: n=%d customers, Q=%g, NV=%d (%.1fs to build distance matrix / road categories)\n', ...
                prob.n, prob.Q, prob.NV, toc(tLoad));
        end

        rc = struct('numRuns',cfg.numRuns, 'referenceAlgo',cfg.referenceAlgo, ...
                    'verbose',cfg.verbose, ...
                    'nsga',struct('popSize',cfg.popSize,'numGen',cfg.numIter,'verbose',false));
        stats = run30_mozoa_paper(prob, rc);

        save(resultFile, 'stats', 'prob', '-v7');
    end

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

    save('large_scale_allresults.mat', 'allResults', '-v7');
    try
        exportLargeScaleSummary(allResults, 'large_scale_allresults.xlsx');
    catch ME
        if cfg.verbose, fprintf('  (Summary Excel export skipped: %s)\n', ME.message); end
    end
end

if cfg.verbose
    fprintf('\n##############################################################\n');
    fprintf('# Sweep complete: %d instance(s) in %.1f min\n', numel(allResults), toc(sweepStart)/60);
    fprintf('# Per-instance files: large_scale_results_<instance>.mat/.xlsx\n');
    fprintf('# Combined summary:   large_scale_allresults.mat/.xlsx\n');
    fprintf('##############################################################\n');
end
end

% =========================================================================
function n = detectGoldenBlockCount(fpath)
% Count how many "[k] with N Customers" instance headers a combined-format
% Golden file contains (kelly-pbs.txt style). Returns 1 if NO such header
% is found, since that means the file is in the single-instance format
% (kelly01.txt..kelly20.txt style) instead, which always contains exactly
% one instance.
fid = fopen(fpath, 'r');
if fid == -1
    error('detectGoldenBlockCount:fileNotFound', 'Could not open "%s".', fpath);
end
raw = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, {'\r\n','\n'});
n = 0;
for i = 1:numel(lines)
    if ~isempty(regexp(lines{i}, '\[(\d+)\]\s*with\s*(\d+)\s*Customers', 'once'))
        n = n + 1;
    end
end
if n == 0
    n = 1;  % single-instance format: exactly one instance, no combined-format header
end
end

% =========================================================================
function exportLargeScaleSummary(allResults, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
if isempty(allResults), return; end
algos = allResults(1).algos;

header = [{'Instance','n'}, algos];
rowsHV = cell(numel(allResults), numel(header));
rowsNDS = cell(numel(allResults), numel(header));
rowsCPU = cell(numel(allResults), numel(header));
for k = 1:numel(allResults)
    r = allResults(k);
    rowsHV(k,:)  = [{r.name, r.n}, num2cell(r.hvMean)];
    rowsNDS(k,:) = [{r.name, r.n}, num2cell(r.ndsMean)];
    rowsCPU(k,:) = [{r.name, r.n}, num2cell(r.cpuMean)];
end
writecell([header; rowsHV],  xlsxFile, 'Sheet', 'HV_mean');
writecell([header; rowsNDS], xlsxFile, 'Sheet', 'NDS_mean');
writecell([header; rowsCPU], xlsxFile, 'Sheet', 'CPU_mean');
end

function o = mergeOptsLarge(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
