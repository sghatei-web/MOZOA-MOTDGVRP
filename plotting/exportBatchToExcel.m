function exportBatchToExcel(batch, xlsxFile)
%EXPORTBATCHTOEXCEL  Write the multi-instance "batch" struct array produced by
%   run_batch_experiments.m to a multi-sheet Excel workbook (.xlsx), as a
%   drop-in companion to (not a replacement for) batch_summary.mat.
%
%   exportBatchToExcel(batch, xlsxFile)
%
%   batch is a struct array, one element per instance, each with fields
%   name, stats, hvMean, hvStd, ndsMean, cpuMean, pWilcoxon (see
%   run_batch_experiments.m). This differs in shape from the per-run
%   algorithm/condition matrices written by exportStatsToExcel.m, so it has
%   its own writer:
%
%   Sheets written:
%     Overview        : one row per instance, the same summary columns
%                        printed to the console by run_batch_experiments.m
%                        (HV/NDS mean per algorithm, Wilcoxon p vs NSGA-2)
%     <InstanceName>   : one sheet per instance, containing that instance's
%                        full stats struct written with the same layout as
%                        exportStatsToExcel.m (Raw_HV/Raw_NDS/Raw_CPU/Summary/
%                        Coverage/Wilcoxon all collapsed into one sheet per
%                        instance, since Excel sheet names must be unique and
%                        short); this is a condensed view, not a full
%                        substitute for calling exportStatsToExcel on that
%                        instance's stats directly if you need the raw
%                        per-run matrices for one particular instance.
%
%   Uses only built-in writecell (no toolbox required; requires MATLAB
%   R2019a or later, when writecell was introduced).
%
%   USAGE
%     batch = run_batch_experiments(cfg);
%     exportBatchToExcel(batch, 'batch_summary.xlsx');

if nargin < 2 || isempty(xlsxFile)
    xlsxFile = 'batch_summary.xlsx';
end
if exist(xlsxFile,'file'), delete(xlsxFile); end

if isempty(batch)
    warning('exportBatchToExcel:empty', 'batch is empty; nothing to export.');
    return;
end

% ---- Overview sheet: one row per instance ----
header = {'Instance','HV_MOZOA','HV_NSGA2','HV_MLNSGA2', ...
          'NDS_MOZOA','NDS_NSGA2','NDS_MLNSGA2', ...
          'CPU_MOZOA','CPU_NSGA2','CPU_MLNSGA2', ...
          'p_MOZOA_vs_NSGA2','p_MOZOA_vs_MLNSGA2'};
rows = cell(numel(batch), numel(header));
for k = 1:numel(batch)
    b = batch(k);
    rows(k,:) = {b.name, ...
        b.hvMean(1), b.hvMean(2), b.hvMean(3), ...
        b.ndsMean(1), b.ndsMean(2), b.ndsMean(3), ...
        b.cpuMean(1), b.cpuMean(2), b.cpuMean(3), ...
        b.pWilcoxon(1), b.pWilcoxon(2)};
end
writecell([header; rows], xlsxFile, 'Sheet', 'Overview');

% ---- one condensed sheet per instance, reusing each instance's full stats -
for k = 1:numel(batch)
    b = batch(k);
    if ~isfield(b,'stats') || isempty(b.stats), continue; end
    sheetName = safeSheetName(b.name, k);
    writeInstanceSheet(xlsxFile, sheetName, b.stats);
end

fprintf('Excel batch results written to %s (%d instance(s)).\n', xlsxFile, numel(batch));
end

% =========================================================================
function writeInstanceSheet(xlsxFile, sheetName, stats)
% Condensed single-sheet dump of one instance's stats: summary table first,
% then coverage/Wilcoxon below it, all in one sheet (Excel has no nested
% sheets, so a full per-instance breakdown uses exportStatsToExcel directly).
rows = {};
rows(end+1,:) = {'--- Summary (best/worst/mean/std/median) ---','','','','',''}; %#ok<AGROW>
rows(end+1,:) = {'Metric','Algorithm','Best','Worst','Mean','Std'}; %#ok<AGROW>
if isfield(stats,'summary')
    labels = {};
    if isfield(stats,'algos'), labels = stats.algos;
    elseif isfield(stats,'conditions'), labels = stats.conditions; end
    metricNames = fieldnames(stats.summary);
    for m = 1:numel(metricNames)
        mName = metricNames{m};
        arr = stats.summary.(mName);
        for a = 1:numel(arr)
            lbl = sprintf('col%d', a);
            if a <= numel(labels), lbl = labels{a}; end
            rows(end+1,:) = {mName, lbl, arr(a).best, arr(a).worst, arr(a).mean, arr(a).std}; %#ok<AGROW>
        end
    end
end
rows(end+1,:) = {'','','','','',''}; %#ok<AGROW>
rows(end+1,:) = {'--- Coverage / Contribution ---','','','','',''}; %#ok<AGROW>
if isfield(stats,'cMetric')
    fn = fieldnames(stats.cMetric);
    for i = 1:numel(fn)
        rows(end+1,:) = {fn{i}, stats.cMetric.(fn{i}), '', '', '', ''}; %#ok<AGROW>
    end
end
if isfield(stats,'contribution')
    fn = fieldnames(stats.contribution);
    for i = 1:numel(fn)
        rows(end+1,:) = {fn{i}, stats.contribution.(fn{i}), '', '', '', ''}; %#ok<AGROW>
    end
end
rows(end+1,:) = {'','','','','',''}; %#ok<AGROW>
rows(end+1,:) = {'--- Wilcoxon p-values ---','','','','',''}; %#ok<AGROW>
if isfield(stats,'wilcoxon')
    fn = fieldnames(stats.wilcoxon);
    for i = 1:numel(fn)
        rows(end+1,:) = {fn{i}, stats.wilcoxon.(fn{i}), '', '', '', ''}; %#ok<AGROW>
    end
end
writecell(rows, xlsxFile, 'Sheet', sheetName);
end

function s = safeSheetName(name, idx)
% Excel sheet names: max 31 chars, no \/?*[]:  -- sanitise and de-duplicate
% by prefixing the instance index if truncation could collide.
s = regexprep(name, '[\\/\?\*\[\]:]', '_');
if numel(s) > 28
    s = s(1:28);
end
s = sprintf('%02d_%s', idx, s);
if numel(s) > 31
    s = s(1:31);
end
end
