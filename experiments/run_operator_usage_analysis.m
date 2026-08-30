function results = run_operator_usage_analysis(cfg)
%RUN_OPERATOR_USAGE_ANALYSIS  Per-operator usage statistics for MOZOA's
%   full (fixed-schedule) configuration, addressing the request for a
%   more direct justification of the fixed operator-selection
%   probabilities than "these operators were more effective".
%
%   results = run_operator_usage_analysis(cfg)
%
%   For each instance and each of MOZOA's eight operators, this reports,
%   pooled over cfg.numRuns independent runs:
%     Frequency          how often the operator was selected, as a
%                         fraction of all operator-selection events
%                         (should closely match the fixed schedule's
%                         assigned probability -- reported mainly as a
%                         sanity check that the RNG draw matches the
%                         intended weights)
%     AcceptanceRate      of the times this operator was applied, the
%                         fraction where the resulting solution was
%                         accepted (replaced the parent), under MOZOA's
%                         two-phase foraging/defence rule
%     ArchiveEntryRate    of the times this operator was applied, the
%                         fraction where the resulting solution actually
%                         entered the external Pareto archive (a
%                         stricter, quality-focused signal than plain
%                         acceptance, since a solution can be "accepted"
%                         under the foraging-phase soft-improvement rule
%                         without being good enough for the archive)
%
%   This uses solveZOA8opAblation.m's 'full' condition (identical to
%   solveZOA8op.m, MOZOA exactly as reported in the paper) which already
%   instruments every operator-selection event; no separate solver is
%   needed.
%
%   INPUTS (all optional; cfg is a struct with any of these fields)
%     cfg.instances   cell array of .vrp filenames (default: same 3
%                     representative instances as run_ablation_study.m)
%     cfg.numRuns     independent runs per instance (default: 30)
%     cfg.popSize     population size (default: 100)
%     cfg.numIter     iteration budget (default: 100)
%     cfg.Amax        archive size (default: 100)
%     cfg.p1          defence probability (default: 0.60)
%     cfg.verbose     default true
%
%   OUTPUT
%     results : struct array, one row per (instance, operator), with
%       fields: instance, n, operator, operatorName, frequency,
%       acceptanceRate, archiveEntryRate
%     Also saved to operator_usage_results.mat/.xlsx.
%
%   USAGE
%     run_operator_usage_analysis;
%     c.numRuns = 30; run_operator_usage_analysis(c);

if nargin < 1 || isempty(cfg), cfg = struct; end
def = struct( ...
    'instances', {{'E-n23-k3.vrp','E-n51-k5.vrp','E-n101-k8.vrp'}}, ...
    'numRuns', 30, 'popSize', 100, 'numIter', 100, 'Amax', 100, 'p1', 0.60, ...
    'verbose', true);
cfg = mergeOptsOpUsage(def, cfg);

opNames = {'Swap (O1)','Insertion (O2)','2-opt (O3)','3-opt (O4)', ...
    'Segment relocation (O5)','OX-archive (O6)','OX-leader (O7)','Double-bridge (O8)'};

results = struct('instance',{},'n',{},'operator',{},'operatorName',{}, ...
    'frequency',{},'acceptanceRate',{},'archiveEntryRate',{});

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

    totalTried    = zeros(1,8);
    totalAccepted = zeros(1,8);
    totalArchived = zeros(1,8);

    for r = 1:cfg.numRuns
        opts = struct('popSize',cfg.popSize, 'numIter',cfg.numIter, ...
            'Amax',cfg.Amax, 'p1',cfg.p1, 'variant','full', ...
            'verbose',false, 'seed',r);
        res = solveZOA8opAblation(prob, opts);
        totalTried    = totalTried    + res.opTried;
        totalAccepted = totalAccepted + res.opAccepted;
        totalArchived = totalArchived + res.opArchived;

        if cfg.verbose && mod(r,10)==0
            fprintf('  run %d/%d done\n', r, cfg.numRuns);
        end
    end

    grandTotal = sum(totalTried);
    for a = 1:8
        freq = totalTried(a) / max(grandTotal, 1);
        accRate = totalAccepted(a) / max(totalTried(a), 1);
        archRate = totalArchived(a) / max(totalTried(a), 1);

        rec = struct('instance',prob.name, 'n',prob.n, 'operator',a, ...
            'operatorName',opNames{a}, 'frequency',freq, ...
            'acceptanceRate',accRate, 'archiveEntryRate',archRate);
        results(end+1) = rec; %#ok<AGROW>

        if cfg.verbose
            fprintf('    %-24s freq=%.3f  accept=%.3f  archiveEntry=%.3f\n', ...
                opNames{a}, freq, accRate, archRate);
        end
    end

    save('operator_usage_results.mat', 'results', '-v7');
    try
        exportOpUsageToExcel(results, 'operator_usage_results.xlsx');
    catch ME
        if cfg.verbose, fprintf('(Excel export skipped: %s)\n', ME.message); end
    end
end

if cfg.verbose
    fprintf('\nOperator usage analysis complete in %.1f min.\n', toc(sweepStart)/60);
    fprintf('Results saved to operator_usage_results.mat/.xlsx\n');
end
end

% =========================================================================
function exportOpUsageToExcel(results, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
header = {'Instance','n','Operator','OperatorName','Frequency', ...
    'AcceptanceRate','ArchiveEntryRate'};
rows = cell(numel(results), numel(header));
for i = 1:numel(results)
    r = results(i);
    rows(i,:) = {r.instance, r.n, r.operator, r.operatorName, ...
        r.frequency, r.acceptanceRate, r.archiveEntryRate};
end
writecell([header; rows], xlsxFile, 'Sheet', 'OperatorUsage');
end

function o = mergeOptsOpUsage(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
