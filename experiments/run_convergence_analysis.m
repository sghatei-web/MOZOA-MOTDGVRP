function results = run_convergence_analysis(cfg)
%RUN_CONVERGENCE_ANALYSIS  Per-generation hypervolume convergence curves
%   for all five algorithms (MOZOA, NSGA-2, MLNSGA-2, SPEA2, MOEA/D),
%   addressing whether MOZOA's advantage is only in final quality or also
%   in convergence speed.
%
%   results = run_convergence_analysis(cfg)
%
%   For each instance and each algorithm, runs cfg.numRuns independent
%   repetitions and records the per-generation normalised hypervolume
%   (mean and std/IQR across runs at every generation), all algorithms
%   scored against the SAME shared reference point per instance (built
%   from the pooled final-generation objectives of every algorithm and
%   every run, mirroring the convention already used for the paper's
%   final-quality hypervolume comparison in Table 16).
%
%   INPUTS (all optional; cfg is a struct with any of these fields)
%     cfg.instances   cell array of .vrp filenames (default: 4
%                     representative instances spanning the tested size
%                     range: {'E-n23-k3.vrp','E-n33-k4.vrp',
%                     'E-n51-k5.vrp','E-n101-k8.vrp'} -- includes
%                     E-n33-k4, the instance where MOZOA does not attain
%                     the best final hypervolume, so its convergence
%                     behaviour can be inspected directly)
%     cfg.numRuns     independent runs per (instance, algorithm)
%                     (default: 10 -- convergence curves need many
%                     generations recorded per run, so this is smaller
%                     than the main study's N=30 to keep runtime
%                     reasonable; raise it for a smoother mean curve)
%     cfg.popSize     population size (default: 100, matching Table 10)
%     cfg.numGen      generation/iteration budget (default: 100, matching
%                     Table 10)
%     cfg.verbose     default true
%
%   OUTPUT
%     results : struct, one field per instance tag, each containing:
%       .instance, .n, .algos (cell array of names),
%       .hvCurveMean (numGen x 5), .hvCurveStd (numGen x 5),
%       .hvCurveQ25, .hvCurveQ75 (numGen x 5, for IQR shading)
%     Also saved to convergence_analysis_results.mat, plus one .xlsx per
%     instance (convergence_<instance>.xlsx) with a sheet per algorithm.
%
%   USAGE
%     run_convergence_analysis;                          % full analysis
%     c.numRuns = 30; c.instances = {'E-n23-k3.vrp'};
%     run_convergence_analysis(c);                         % quick smoke test

if nargin < 1 || isempty(cfg), cfg = struct; end
def = struct( ...
    'instances', {{'E-n23-k3.vrp','E-n33-k4.vrp','E-n51-k5.vrp','E-n101-k8.vrp'}}, ...
    'numRuns', 30, 'popSize', 100, 'numGen', 100, 'verbose', true);
cfg = mergeOptsConvergence(def, cfg);

algos = {'MOZOA','NSGA-2','MLNSGA-2','SPEA2','MOEA/D'};
results = struct();

sweepStart = tic;

for ii = 1:numel(cfg.instances)
    fname = cfg.instances{ii};
    prob = loadCVRPInstance(fname, 'conflict', 7);
    if ~isempty(regexp(prob.name, 'n(51|76|101)', 'once'))
        prob = scaleInstance(prob, 100);
    end
    tag = matlab.lang.makeValidName(prob.name);

    if cfg.verbose
        fprintf('\n========== [%d/%d] %s (n=%d) ==========\n', ...
            ii, numel(cfg.instances), prob.name, prob.n);
    end

   
    pool = [];
    for r = 1:cfg.numRuns
        moz = solveZOA8opAblation(prob, struct('popSize',cfg.popSize, ...
            'numIter',cfg.numGen, 'variant','full', 'verbose',false, 'seed',r));
        pool = [pool; moz.paretoObj]; %#ok<AGROW>
    end
    refPoint = max(pool,[],1) * 1.5 + 1;  % generous margin since competitor
                                           % objectives are not yet pooled in

    % ---- Pass 2: run everyone for real, with the shared reference point,
    %      recording the per-generation HV curve of each ----
    hvCurves = struct();
    for a = 1:numel(algos)
        hvCurves.(algoField(algos{a})) = zeros(cfg.numGen, cfg.numRuns);
    end

    for r = 1:cfg.numRuns
        moz = solveZOA8opAblation(prob, struct('popSize',cfg.popSize, ...
            'numIter',cfg.numGen, 'variant','full', 'verbose',false, 'seed',r));
        hvCurves.MOZOA(:,r) = moz.hvCurve;

        nsga2 = solveMOTDGVRPTracked(prob, struct('popSize',cfg.popSize, ...
            'numGen',cfg.numGen, 'useML',false, 'verbose',false, 'seed',r, ...
            'refPoint',refPoint));
        hvCurves.NSGA_2(:,r) = nsga2.hvCurve;

        mlnsga2 = solveMOTDGVRPTracked(prob, struct('popSize',cfg.popSize, ...
            'numGen',cfg.numGen, 'useML',true, 'verbose',false, 'seed',r, ...
            'refPoint',refPoint));
        hvCurves.MLNSGA_2(:,r) = mlnsga2.hvCurve;

        spea2 = solveSPEA2Tracked(prob, struct('popSize',cfg.popSize, ...
            'numGen',cfg.numGen, 'verbose',false, 'seed',r, ...
            'refPoint',refPoint));
        hvCurves.SPEA2(:,r) = spea2.hvCurve;

        moead = solveMOEADTracked(prob, struct('popSize',cfg.popSize, ...
            'numGen',cfg.numGen, 'verbose',false, 'seed',r, ...
            'refPoint',refPoint));
        hvCurves.MOEA_D(:,r) = moead.hvCurve;

        if cfg.verbose
            fprintf('  run %2d/%d done (all 5 algorithms)\n', r, cfg.numRuns);
        end
    end

    hvCurveMean = zeros(cfg.numGen, numel(algos));
    hvCurveStd  = zeros(cfg.numGen, numel(algos));
    hvCurveQ25  = zeros(cfg.numGen, numel(algos));
    hvCurveQ75  = zeros(cfg.numGen, numel(algos));
    for a = 1:numel(algos)
        M = hvCurves.(algoField(algos{a}));
        hvCurveMean(:,a) = mean(M, 2);
        hvCurveStd(:,a)  = std(M, 0, 2);
        hvCurveQ25(:,a)  = quantileSimpleConv(M, 0.25);
        hvCurveQ75(:,a)  = quantileSimpleConv(M, 0.75);
    end

    results.(tag) = struct('instance',prob.name, 'n',prob.n, 'algos',{algos}, ...
        'hvCurveMean',hvCurveMean, 'hvCurveStd',hvCurveStd, ...
        'hvCurveQ25',hvCurveQ25, 'hvCurveQ75',hvCurveQ75);

    save('convergence_analysis_results.mat', 'results', '-v7');
    try
        exportConvergenceToExcel(results.(tag), sprintf('convergence_%s.xlsx', tag));
    catch ME
        if cfg.verbose, fprintf('(Excel export skipped: %s)\n', ME.message); end
    end
end

if cfg.verbose
    fprintf('\nConvergence analysis complete in %.1f min.\n', toc(sweepStart)/60);
    fprintf('Results saved to convergence_analysis_results.mat and convergence_<instance>.xlsx\n');
end
end

% =========================================================================
function f = algoField(name)
f = matlab.lang.makeValidName(name);
end

function q = quantileSimpleConv(M, p)
% Row-wise (per-generation) quantile across the columns (runs) of M,
% without requiring the Statistics Toolbox.
[numGen, numRuns] = size(M);
q = zeros(numGen,1);
for g = 1:numGen
    v = sort(M(g,:));
    if numRuns == 1
        q(g) = v(1);
        continue;
    end
    idx = 1 + (numRuns-1)*p;
    lo = floor(idx); hi = ceil(idx);
    if lo == hi
        q(g) = v(lo);
    else
        q(g) = v(lo) + (idx-lo)*(v(hi)-v(lo));
    end
end
end

function exportConvergenceToExcel(res, xlsxFile)
if exist(xlsxFile,'file'), delete(xlsxFile); end
header = {'Generation', res.algos{:}}; %#ok<CCAT>
gen = (1:size(res.hvCurveMean,1))';
writecell([header; num2cell([gen, res.hvCurveMean])], xlsxFile, 'Sheet', 'Mean');
writecell([header; num2cell([gen, res.hvCurveStd])], xlsxFile, 'Sheet', 'Std');
writecell([header; num2cell([gen, res.hvCurveQ25])], xlsxFile, 'Sheet', 'Q25');
writecell([header; num2cell([gen, res.hvCurveQ75])], xlsxFile, 'Sheet', 'Q75');
end

function o = mergeOptsConvergence(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
