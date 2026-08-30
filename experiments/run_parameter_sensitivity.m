function results = run_parameter_sensitivity(cfg)
%RUN_PARAMETER_SENSITIVITY  Parameter sensitivity analysis for MOZOA
%   (Reviewer request: demonstrate robustness to NP, Tmax, and p1).
%
%   results = run_parameter_sensitivity(cfg)
%
%   For each of three representative instances (small/medium/large by
%   default: E-n23-k3, E-n51-k5, E-n101-k8), this varies ONE parameter at
%   a time around MOZOA's paper defaults (NP=100, Tmax=100, p1=0.60)
%   while holding the other two fixed, runs cfg.numRuns independent
%   repetitions per setting, and records mean/std hypervolume and mean
%   CPU time. This directly answers the reviewer's request for a
%   parameter sensitivity study rather than a single fixed configuration.
%
%   INPUTS (all optional; cfg is a struct with any of these fields)
%     cfg.instances   cell array of .vrp filenames to test
%                     (default: {'E-n23-k3.vrp','E-n51-k5.vrp','E-n101-k8.vrp'})
%     cfg.numRuns     independent runs per (instance, parameter, value)
%                     combination (default: 10 -- raise to match your
%                     main experiments' run count if you want directly
%                     comparable precision)
%     cfg.npValues    population-size values to test, holding
%                     Tmax=100, p1=0.60 fixed (default: [50 100 150 200])
%     cfg.tmaxValues  iteration-budget values to test, holding
%                     NP=100, p1=0.60 fixed (default: [50 100 150 200])
%     cfg.p1Values    defence-probability values to test, holding
%                     NP=100, Tmax=100 fixed (default: [0.3 0.45 0.6 0.75 0.9])
%     cfg.verbose     default true
%
%   OUTPUT
%     results : struct with fields .NP, .Tmax, .p1, each a struct array
%       with one row per (instance, value) combination:
%         .instance, .value, .hvMean, .hvStd, .cpuMean, .cpuStd
%     Also saved to parameter_sensitivity_results.mat/.xlsx.
%
%   USAGE
%     run_parameter_sensitivity;                       % full default sweep
%     c.numRuns = 30; run_parameter_sensitivity(c);
%     c.instances = {'E-n23-k3.vrp'}; c.numRuns = 30;
%     run_parameter_sensitivity(c);                     % single instance, fast

if nargin < 1 || isempty(cfg), cfg = struct; end
def = struct( ...
    'instances', {{'E-n23-k3.vrp','E-n51-k5.vrp','E-n101-k8.vrp'}}, ...
    'numRuns', 30, ...
    'npValues', [50 100 150 200], ...
    'tmaxValues', [50 100 150 200], ...
    'p1Values', [0.3 0.45 0.6 0.75 0.9], ...
    'verbose', true);
cfg = mergeOptsSens(def, cfg);

% base (paper-default) settings that stay fixed unless the swept
% parameter is itself the one being varied
BASE_NP   = 100;
BASE_TMAX = 100;
BASE_P1   = 0.60;
BASE_AMAX = 100;

results.NP   = struct('instance',{},'value',{},'hvMean',{},'hvStd',{},'cpuMean',{},'cpuStd',{});
results.Tmax = struct('instance',{},'value',{},'hvMean',{},'hvStd',{},'cpuMean',{},'cpuStd',{});
results.p1   = struct('instance',{},'value',{},'hvMean',{},'hvStd',{},'cpuMean',{},'cpuStd',{});

sweepStart = tic;

for i = 1:numel(cfg.instances)
    fname = cfg.instances{i};
    prob = loadCVRPInstance(fname, 'conflict', 7);
    if ~isempty(regexp(prob.name, 'n(51|76|101)', 'once'))
        prob = scaleInstance(prob, 100);
    end

    if cfg.verbose
        fprintf('\n========== Instance: %s (n=%d) ==========\n', prob.name, prob.n);
    end

    % ---------- sweep NP (population size) ----------
    for np = cfg.npValues
        [hvMean, hvStd, cpuMean, cpuStd] = runSweepPoint(prob, np, BASE_TMAX, BASE_AMAX, BASE_P1, cfg.numRuns, cfg.verbose, 'NP', np);
        rec = struct('instance',prob.name,'value',np,'hvMean',hvMean,'hvStd',hvStd,'cpuMean',cpuMean,'cpuStd',cpuStd);
        results.NP(end+1) = rec; %#ok<AGROW>
    end

    % ---------- sweep Tmax (iteration budget) ----------
    for tmax = cfg.tmaxValues
        [hvMean, hvStd, cpuMean, cpuStd] = runSweepPoint(prob, BASE_NP, tmax, BASE_AMAX, BASE_P1, cfg.numRuns, cfg.verbose, 'Tmax', tmax);
        rec = struct('instance',prob.name,'value',tmax,'hvMean',hvMean,'hvStd',hvStd,'cpuMean',cpuMean,'cpuStd',cpuStd);
        results.Tmax(end+1) = rec; %#ok<AGROW>
    end

    % ---------- sweep p1 (defence probability) ----------
    for p1 = cfg.p1Values
        [hvMean, hvStd, cpuMean, cpuStd] = runSweepPoint(prob, BASE_NP, BASE_TMAX, BASE_AMAX, p1, cfg.numRuns, cfg.verbose, 'p1', p1);
        rec = struct('instance',prob.name,'value',p1,'hvMean',hvMean,'hvStd',hvStd,'cpuMean',cpuMean,'cpuStd',cpuStd);
        results.p1(end+1) = rec; %#ok<AGROW>
    end

    % checkpoint after every instance so a long sweep is resumable in spirit
    % (re-running redoes everything, but results so far are saved to disk)
    save('parameter_sensitivity_results.mat', 'results', '-v7');
    try
        exportSensitivityToExcel(results, 'parameter_sensitivity_results.xlsx');
    catch ME
        if cfg.verbose, fprintf('(Excel export skipped: %s)\n', ME.message); end
    end
end

if cfg.verbose
    fprintf('\nParameter sensitivity sweep complete in %.1f min.\n', toc(sweepStart)/60);
    fprintf('Results saved to parameter_sensitivity_results.mat/.xlsx\n');
end
end

% =========================================================================
function [hvMean, hvStd, cpuMean, cpuStd] = runSweepPoint(prob, np, tmax, amax, p1, numRuns, verbose, paramName, paramValue)
hv = zeros(numRuns,1);
cpu = zeros(numRuns,1);
opts = struct('popSize',np,'numIter',tmax,'Amax',amax,'p1',p1,'verbose',false);
for r = 1:numRuns
    opts.seed = r;
    res = solveZOA8op(prob, opts);
    ref = max(res.paretoObj,[],1)*1.1 + 1;
    hv(r) = hypervolume(res.paretoObj, ref);
    cpu(r) = res.cpuTime;
end
hvMean = mean(hv); hvStd = std(hv);
cpuMean = mean(cpu); cpuStd = std(cpu);
if verbose
    fprintf('  %-6s = %-6g : HV mean=%.4g (std=%.4g), CPU mean=%.2fs (std=%.2fs)\n', ...
        paramName, paramValue, hvMean, hvStd, cpuMean, cpuStd);
end
end

function exportSensitivityToExcel(results, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
fields = {'NP','Tmax','p1'};
for f = 1:numel(fields)
    fn = fields{f};
    R = results.(fn);
    if isempty(R), continue; end
    header = {'Instance','Value','HV_mean','HV_std','CPU_mean','CPU_std'};
    rows = cell(numel(R), numel(header));
    for i = 1:numel(R)
        rows(i,:) = {R(i).instance, R(i).value, R(i).hvMean, R(i).hvStd, R(i).cpuMean, R(i).cpuStd};
    end
    writecell([header; rows], xlsxFile, 'Sheet', fn);
end
end

function o = mergeOptsSens(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
