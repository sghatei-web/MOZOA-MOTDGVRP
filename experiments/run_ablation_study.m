function results = run_ablation_study(cfg)
%RUN_ABLATION_STUDY  Ablation study isolating the contribution of MOZOA's
%   external Pareto archive, recombination operators (O6/O7), and fixed
%   operator-selection schedule.
%
%   results = run_ablation_study(cfg)
%
%   Runs four conditions, all otherwise identical (same population
%   initialisation, same 2-opt/3-opt/swap/insertion/segment-relocation/
%   double-bridge operators, same acceptance rule where applicable):
%
%     full             MOZOA exactly as reported in the paper (archive +
%                       all 8 operators + fixed schedule)
%     noArchive        the external Pareto archive is removed; only the
%                       current population's own non-dominated subset is
%                       available for leader/partner selection and HV
%                       tracking each iteration (no memory across
%                       iterations beyond the population)
%     noRecombination  operators O6 (archive OX) and O7 (leader OX) are
%                       removed from the operator set; the remaining six
%                       operators are drawn with the fixed schedule's
%                       relative weights renormalised to sum to 1
%     uniformSelection archive + all 8 operators kept, but operators are
%                       drawn uniformly at random each step instead of
%                       the fixed non-uniform schedule
%
%   For each condition, on each instance, cfg.numRuns independent runs
%   are performed and HV, NDS, IGD+, the C-metric of that condition
%   against 'full', and the Contribution of that condition against
%   'full' are reported (IGD+ uses a reference set pooled across all four
%   conditions' fronts on that instance -- see buildReferenceSet.m).
%
%   INPUTS (all optional; cfg is a struct with any of these fields)
%     cfg.instances   cell array of .vrp filenames to test
%                     (default: {'E-n23-k3.vrp','E-n51-k5.vrp','E-n101-k8.vrp'},
%                     matching the parameter-sensitivity study's
%                     representative small/medium/large instances)
%     cfg.numRuns     independent runs per (instance, condition)
%                     (default: 30, matching the paper's main sample size)
%     cfg.popSize     population size (default: 100, matching Table 10)
%     cfg.numIter     iteration budget (default: 100, matching Table 10)
%     cfg.Amax        archive size (default: 100, matching Table 10)
%     cfg.p1          defence probability (default: 0.60, matching Table 10)
%     cfg.verbose     default true
%
%   OUTPUT
%     results : struct array, one row per (instance, condition), with
%       fields: instance, n, condition, hvMean, hvStd, ndsMean, ndsStd,
%       igdMean, igdStd, cMetricVsFull, cMetricFullVsThis,
%       contributionVsFull, cpuMean
%     Also saved to ablation_study_results.mat/.xlsx.
%
%   USAGE
%     run_ablation_study;                              % full study
%     c.numRuns = 30; c.instances = {'E-n23-k3.vrp'};
%     run_ablation_study(c);

if nargin < 1 || isempty(cfg), cfg = struct; end
def = struct( ...
    'instances', {{'E-n23-k3.vrp','E-n51-k5.vrp','E-n101-k8.vrp'}}, ...
    'numRuns', 30, 'popSize', 100, 'numIter', 100, 'Amax', 100, 'p1', 0.60, ...
    'verbose', true);
cfg = mergeOptsAblationStudy(def, cfg);

conditions = {'full','noArchive','noRecombination','uniformSelection'};
condLabels = struct('full','MOZOA-full', 'noArchive','MOZOA-noArchive', ...
    'noRecombination','MOZOA-noRecomb', 'uniformSelection','MOZOA-uniformSel');

results = struct('instance',{},'n',{},'condition',{},'hvMean',{},'hvStd',{}, ...
    'ndsMean',{},'ndsStd',{},'igdMean',{},'igdStd',{}, ...
    'cMetricVsFull',{},'cMetricFullVsThis',{},'contributionVsFull',{},'cpuMean',{});

sweepStart = tic;

for ii = 1:numel(cfg.instances)
    fname = cfg.instances{ii};
    prob = loadCVRPInstance(fname, 'conflict', 7);
    if ~isempty(regexp(prob.name, 'n(51|76|101)', 'once'))
        prob = scaleInstance(prob, 100);
    end

    if cfg.verbose
        fprintf('\n========== [%d/%d] %s (n=%d) ==========\n', ...
            ii, numel(cfg.instances), prob.name, prob.n);
    end

    fronts = struct();
    hvAll  = struct();
    ndsAll = struct();
    cpuAll = struct();

    for c = 1:numel(conditions)
        cond = conditions{c};
        hvAll.(cond)  = zeros(cfg.numRuns,1);
        ndsAll.(cond) = zeros(cfg.numRuns,1);
        cpuAll.(cond) = zeros(cfg.numRuns,1);
        fronts.(cond) = cell(cfg.numRuns,1);

        for r = 1:cfg.numRuns
            opts = struct('popSize',cfg.popSize, 'numIter',cfg.numIter, ...
                'Amax',cfg.Amax, 'p1',cfg.p1, 'variant',cond, ...
                'verbose',false, 'seed',r);
            res = solveZOA8opAblation(prob, opts);
            fronts.(cond){r} = res.paretoObj;
            ndsAll.(cond)(r) = size(res.paretoObj,1);
            cpuAll.(cond)(r) = res.cpuTime;
        end

        if cfg.verbose
            fprintf('  %-20s done (%d runs)\n', condLabels.(cond), cfg.numRuns);
        end
    end

    pool = [];
    for c = 1:numel(conditions)
        pool = [pool; vertcat(fronts.(conditions{c}){:})]; %#ok<AGROW>
    end
    refPoint = max(pool,[],1)*1.1 + 1;

    for c = 1:numel(conditions)
        cond = conditions{c};
        for r = 1:cfg.numRuns
            hvAll.(cond)(r) = hypervolume(fronts.(cond){r}, refPoint);
        end
    end

    allFrontsForRef = {};
    for c = 1:numel(conditions)
        allFrontsForRef = [allFrontsForRef, fronts.(conditions{c})]; %#ok<AGROW>
    end
    R = buildReferenceSet(allFrontsForRef(:)');

    igdAll = struct();
    for c = 1:numel(conditions)
        cond = conditions{c};
        igdAll.(cond) = zeros(cfg.numRuns,1);
        for r = 1:cfg.numRuns
            igdAll.(cond)(r) = igdPlus(fronts.(cond){r}, R);
        end
    end

    [~, bestIdxFull] = max(hvAll.full);
    fullFront = fronts.full{bestIdxFull};

    for c = 1:numel(conditions)
        cond = conditions{c};
        [~, bestIdxCond] = max(hvAll.(cond));
        condFront = fronts.(cond){bestIdxCond};

        cThisVsFull = cMetric(condFront, fullFront);
        cFullVsThis = cMetric(fullFront, condFront);
        contribThisVsFull = contributionMetric(condFront, fullFront);

        rec = struct('instance',prob.name, 'n',prob.n, ...
            'condition',condLabels.(cond), ...
            'hvMean',mean(hvAll.(cond)), 'hvStd',std(hvAll.(cond)), ...
            'ndsMean',mean(ndsAll.(cond)), 'ndsStd',std(ndsAll.(cond)), ...
            'igdMean',mean(igdAll.(cond)), 'igdStd',std(igdAll.(cond)), ...
            'cMetricVsFull',cThisVsFull, 'cMetricFullVsThis',cFullVsThis, ...
            'contributionVsFull',contribThisVsFull, ...
            'cpuMean',mean(cpuAll.(cond)));
        results(end+1) = rec; %#ok<AGROW>

        if cfg.verbose
            fprintf('    %-20s HV=%.4g(%.3g)  NDS=%.2f  IGD+=%.4g  Contrib(vs full)=%.3f\n', ...
                condLabels.(cond), rec.hvMean, rec.hvStd, rec.ndsMean, rec.igdMean, ...
                rec.contributionVsFull);
        end
    end

    save('ablation_study_results.mat', 'results', '-v7');
    try
        exportAblationToExcel(results, 'ablation_study_results.xlsx');
    catch ME
        if cfg.verbose, fprintf('(Excel export skipped: %s)\n', ME.message); end
    end
end

if cfg.verbose
    fprintf('\nAblation study complete in %.1f min.\n', toc(sweepStart)/60);
    fprintf('Results saved to ablation_study_results.mat/.xlsx\n');
end
end

% =========================================================================
function exportAblationToExcel(results, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
header = {'Instance','n','Condition','HV_mean','HV_std','NDS_mean','NDS_std', ...
    'IGDplus_mean','IGDplus_std','C_cond_vs_full','C_full_vs_cond', ...
    'Contribution_vs_full','CPU_mean'};
rows = cell(numel(results), numel(header));
for i = 1:numel(results)
    r = results(i);
    rows(i,:) = {r.instance, r.n, r.condition, r.hvMean, r.hvStd, ...
        r.ndsMean, r.ndsStd, r.igdMean, r.igdStd, r.cMetricVsFull, ...
        r.cMetricFullVsThis, r.contributionVsFull, r.cpuMean};
end
writecell([header; rows], xlsxFile, 'Sheet', 'Ablation');
end

function o = mergeOptsAblationStudy(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
