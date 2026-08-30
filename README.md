# MOZOA for the Multi-Objective Time-Dependent Green Vehicle Routing Problem

This repository contains the MATLAB implementation accompanying the paper:

> **A Multi-Operator Zebra Optimization Algorithm with Pareto Archiving for the Multi-Objective Time-Dependent Green Vehicle Routing Problem**

Repository: https://github.com/sghatei-web/MOZOA-MOTDGVRP

MOZOA solves the MOTDGVRP with three minimization objectives: total distance,
transportation time, and fuel consumption. The repository also contains the
four comparison algorithms used in the paper: NSGA-2, MLNSGA-2, SPEA2, and
MOEA/D.

## Requirements

- MATLAB R2021b or later
- Base MATLAB for the optimization algorithms
- Statistics and Machine Learning Toolbox only for MATLAB's `signrank`
  function used in the Wilcoxon signed-rank analysis
- Python plotting is optional; dependencies are listed in
  `plotting/requirements_plotting.txt`

The experiments reported in the paper use 30 independent runs for every
algorithm-instance or experimental-condition comparison.

## Repository structure

- `src/`: MOZOA, comparison algorithms, objective evaluation, decoding,
  Pareto archiving, and performance indicators
- `experiments/`: scripts for the main comparison and supplementary studies
- `benchmarks/set_e/`: local destination for the Christofides Set E instances
  downloaded from the official benchmark collection (not redistributed here)
- `plotting/`: MATLAB and optional Python plotting/export scripts

## Quick start

Start MATLAB in the repository root and run:

```matlab
addpath(genpath('src'));
addpath(genpath('experiments'));
addpath(genpath('plotting'));

cfg = struct('dataset','E-n23-k3.vrp', ...
             'numRuns',30, ...
             'resume',true, ...
             'verbose',true);
stats = run30_mozoa_paper([], cfg);
```

The default protocol uses population size 100, 100 iterations/generations,
archive size 100 for MOZOA, and paired random seeds across algorithms.

## Primary Set E experiment

Download the 13 Christofides Set E `.vrp` files from the public
CVRPLIB/NEO benchmark collection and place them in `benchmarks/set_e/` before
running the primary experiment. The required filenames are:

```text
E-n13-k4.vrp, E-n22-k4.vrp, E-n23-k3.vrp, E-n30-k3.vrp,
E-n31-k7.vrp, E-n33-k4.vrp, E-n51-k5.vrp, E-n76-k7.vrp,
E-n76-k8.vrp, E-n76-k10.vrp, E-n76-k14.vrp,
E-n101-k8.vrp, E-n101-k14.vrp
```

To run all 13 Set E instances with 30 independent runs:

```matlab
addpath(genpath('src'));
addpath(genpath('experiments'));
addpath(genpath('plotting'));

cfg = struct('numRuns',30, ...
             'referenceAlgo','MOZOA', ...
             'scale100',true, ...
             'verbose',true);
allResults = run_all_instances_for_paper(cfg);
```

The script applies the demand/capacity scaling described in the paper to the
E-n51, E-n76, and E-n101 groups. Checkpoint files allow interrupted runs to be
resumed.

## Supplementary experiments

All supplied experiment scripts default to 30 independent runs:

```matlab
run_ablation_study;
run_convergence_analysis;
run_parameter_sensitivity;
run_operator_usage_analysis;
run_fixed_vs_random_schedule;
```

The large-scale Golden benchmark experiment is started with:

```matlab
cfg = struct('instancesDir','PATH_TO_GOLDEN_INSTANCES', ...
             'numRuns',30, ...
             'popSize',100, ...
             'numIter',100);
results = run_large_scale_evaluation(cfg);
```

Golden benchmark files are not redistributed in this repository. They can be
obtained from the public CVRPLIB/NEO benchmark collection and placed in the
directory passed through `instancesDir`.

## Main implementation entry points

- `src/solveZOA8op.m`: proposed MOZOA solver
- `src/applyOperator.m`: eight discrete search operators
- `src/archiveUpdate.m`: external Pareto archive update
- `src/crowdingDistance.m`: crowding-distance calculation
- `src/evaluate.m`: three-objective MOTDGVRP evaluation
- `src/arcTimeFuel.m`: FIFO travel-time and fuel calculation
- `src/solveMOTDGVRP.m`: NSGA-2 and MLNSGA-2
- `src/solveSPEA2.m`: SPEA2
- `src/solveMOEAD.m`: MOEA/D

## Reproducibility notes

- The optimization methods use paired seeds: run `r` uses seed `r` for each
  compared method.
- Cross-algorithm hypervolume uses a common instance-specific reference point
  computed from pooled fronts.
- The true Pareto fronts are unknown; IGD+ is used only in the closed ablation
  comparison with a pooled reference set.
- CPU times depend on hardware and MATLAB configuration and should be compared
  only under a common platform.

## Citation

Please cite the accompanying paper and the archived Zenodo release:

> Ghatei, S., Kusetogullari, H., Arasteh, B., & TaghipourEivazi, S. (2026).
> *MOZOA for the Multi-Objective Time-Dependent Green Vehicle Routing Problem*
> (Version 1.0.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.22170033

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22170033.svg)](https://doi.org/10.5281/zenodo.22170033)

## License

The source code is released under the MIT License. See `LICENSE`.
