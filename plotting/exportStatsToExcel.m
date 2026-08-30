function exportStatsToExcel(stats, xlsxFile)
%EXPORTSTATSTOEXCEL  Write a MOZOA project "stats" struct to a multi-sheet
%   Excel workbook (.xlsx), as a drop-in companion to (not a replacement for)
%   the .mat file produced by run30_MOZOA_GVRP.m, run30_extended_competitors.m,
%   and run_ablation_study.m.
%
%   exportStatsToExcel(stats, xlsxFile)
%
%   Works generically with the "stats" struct produced by any of the three
%   scripts above: it reads whichever of stats.algos / stats.conditions is
%   present as the column labels, and writes whichever of
%   stats.cMetric / stats.contribution / stats.wilcoxon are present (each is
%   optional, since run_ablation_study.m has no cMetric/contribution).
%
%   Sheets written:
%     Raw_HV     : one row per run, one column per algorithm/condition
%     Raw_NDS    : same layout, number of non-dominated solutions
%     Raw_CPU    : same layout, CPU time (seconds)
%     Summary    : best/worst/mean/std/median of HV, NDS, CPU per algorithm
%     Coverage   : mean C-metric and Contribution vs. MOZOA/full (if present)
%     Wilcoxon   : p-values vs. MOZOA/full (if present)
%
%   Uses only built-in writecell/writematrix (no toolbox required; requires
%   MATLAB R2019a or later, when writecell was introduced), so it runs on any
%   sufficiently recent MATLAB installation without extra dependencies.
%
%   USAGE
%     stats = run30_MOZOA_GVRP(prob);
%     exportStatsToExcel(stats, 'run30_results_E-n23-k3.xlsx');
%
%     stats = run30_extended_competitors(prob);
%     exportStatsToExcel(stats, 'run30ext_results_E-n23-k3.xlsx');
%
%     stats = run_ablation_study(prob);
%     exportStatsToExcel(stats, 'ablation_results_E-n23-k3.xlsx');

if nargin < 2 || isempty(xlsxFile)
    xlsxFile = 'stats_export.xlsx';
end
if exist(xlsxFile,'file'), delete(xlsxFile); end   % avoid mixing with a stale file

% ---- column labels: algos (comparison scripts) or conditions (ablation) ----
if isfield(stats,'algos')
    labels = stats.algos;
elseif isfield(stats,'conditions')
    labels = stats.conditions;
else
    error('exportStatsToExcel:noLabels', ...
        'stats has neither .algos nor .conditions; not a recognised stats struct.');
end
labels = labels(:)';   % row cell array

% ---- raw per-run matrices ----
writeRawSheet(xlsxFile, 'Raw_HV',  stats.HV,  labels);
writeRawSheet(xlsxFile, 'Raw_NDS', stats.NDS, labels);
writeRawSheet(xlsxFile, 'Raw_CPU', stats.CPU, labels);

% ---- summary (best/worst/mean/std/median) ----
if isfield(stats,'summary')
    header = {'Metric','Algorithm','Best','Worst','Mean','Std','Median'};
    rows = {};
    metricNames = fieldnames(stats.summary);
    for m = 1:numel(metricNames)
        mName = metricNames{m};
        arr = stats.summary.(mName);
        for a = 1:numel(arr)
            lbl = safeLabel(labels, a);
            rows(end+1,:) = {mName, lbl, arr(a).best, arr(a).worst, ...
                              arr(a).mean, arr(a).std, arr(a).median}; %#ok<AGROW>
        end
    end
    writecell([header; rows], xlsxFile, 'Sheet', 'Summary');
end

% ---- coverage: C-metric + Contribution (only present for the 3/5-algorithm
%      comparison scripts, not for run_ablation_study) ----
hasCov = isfield(stats,'cMetric') && ~isempty(fieldnames(stats.cMetric));
hasContrib = isfield(stats,'contribution') && ~isempty(fieldnames(stats.contribution));
if hasCov || hasContrib
    header = {'Quantity','Value'};
    rows = {};
    if hasCov
        fn = fieldnames(stats.cMetric);
        for i = 1:numel(fn)
            rows(end+1,:) = {fn{i}, stats.cMetric.(fn{i})}; %#ok<AGROW>
        end
    end
    if hasContrib
        fn = fieldnames(stats.contribution);
        for i = 1:numel(fn)
            rows(end+1,:) = {fn{i}, stats.contribution.(fn{i})}; %#ok<AGROW>
        end
    end
    writecell([header; rows], xlsxFile, 'Sheet', 'Coverage');
end

% ---- Wilcoxon p-values ----
if isfield(stats,'wilcoxon') && ~isempty(fieldnames(stats.wilcoxon))
    header = {'Comparison','p_value','Significant_at_0.05'};
    fn = fieldnames(stats.wilcoxon);
    rows = cell(numel(fn),3);
    for i = 1:numel(fn)
        p = stats.wilcoxon.(fn{i});
        rows(i,:) = {fn{i}, p, ~isnan(p) && p < 0.05};
    end
    writecell([header; rows], xlsxFile, 'Sheet', 'Wilcoxon');
end

% ---- reference point used for HV, for traceability ----
if isfield(stats,'ref')
    writecell([{'Objective_1_dist','Objective_2_time','Objective_3_fuel'}; ...
               num2cell(stats.ref)], xlsxFile, 'Sheet', 'HV_Reference_Point');
end

% ---- best-run Pareto front per algorithm (distance/time/fuel points), for
%      3-D Pareto-front plotting; one sheet per algorithm since front sizes
%      differ across algorithms and cannot share a single rectangular sheet --
%      each sheet is named "PF_<algorithm>" (sanitised to a valid, <=31-char
%      Excel sheet name) ----
if isfield(stats,'bestFront') && ~isempty(stats.bestFront)
    for a = 1:numel(stats.bestFront)
        F = stats.bestFront{a};
        if isempty(F), continue; end
        lbl = safeLabel(labels, a);
        sheetName = ['PF_' safeSheetTag(lbl)];
        header = {'Distance','Time','Fuel'};
        writecell([header; num2cell(F)], xlsxFile, 'Sheet', sheetName);
    end
end

fprintf('Excel results written to %s\n', xlsxFile);
end

% =========================================================================
function writeRawSheet(xlsxFile, sheetName, M, labels)
if isempty(M), return; end
nRuns = size(M,1);
header = [{'Run'}, labels];
runCol = num2cell((1:nRuns)');
body = [runCol, num2cell(M)];
writecell([header; body], xlsxFile, 'Sheet', sheetName);
end

function lbl = safeLabel(labels, idx)
if idx <= numel(labels)
    lbl = labels{idx};
else
    lbl = sprintf('col%d', idx);
end
end

function s = safeSheetTag(name)
% Excel sheet names: max 31 chars total, no \ / ? * [ ] :
s = regexprep(name, '[\\/\?\*\[\]:]', '');
if numel(s) > 27, s = s(1:27); end  % leave room for the "PF_" prefix
end
