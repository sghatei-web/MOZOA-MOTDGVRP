# Golden/Kelly large-scale CVRP benchmark

The large-scale experiments use the 20 Golden--Wasil--Kelly--Chao CVRP
instances, commonly distributed as `kelly01.txt` through `kelly20.txt`.

The benchmark files are **not redistributed in this repository** because the
downloaded archive does not include an explicit redistribution licence. Obtain
the files from the NEO Research Group's public CVRP benchmark page:

- Benchmark catalogue:
  https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-instances/
- Direct archive supplied by that catalogue:
  https://neo.lcc.uma.es/vrp/wp-content/data/instances/kelly/kelly.zip

## Required files

Extract the archive and copy the following files into this directory:

```text
kelly01.txt  kelly02.txt  kelly03.txt  kelly04.txt  kelly05.txt
kelly06.txt  kelly07.txt  kelly08.txt  kelly09.txt  kelly10.txt
kelly11.txt  kelly12.txt  kelly13.txt  kelly14.txt  kelly15.txt
kelly16.txt  kelly17.txt  kelly18.txt  kelly19.txt  kelly20.txt
```

The official archive may also contain:

- `kelly-pbs.txt`: all 20 instances in one combined file;
- `kelly-solns.txt`: reference solutions.

Use **either** the 20 individual files **or** `kelly-pbs.txt`, not both in the
same experiment directory. Otherwise, `run_large_scale_evaluation.m` will
detect the same instances twice. The `kelly-solns.txt` file is not required
to execute the optimization algorithms.

## Run the experiment

From the repository root:

```matlab
addpath(genpath('src'));
addpath(genpath('experiments'));

cfg = struct('instancesDir', fullfile('benchmarks','golden'), ...
             'numRuns', 30, ...
             'popSize', 100, ...
             'numIter', 100);

results = run_large_scale_evaluation(cfg);
```

The supplied parsers support both the individual-file format and the combined
`kelly-pbs.txt` format.

## Benchmark attribution

These instances are attributed to Golden, Wasil, Kelly, and Chao and are
distributed by the NEO Research Group as a 20-instance large-scale CVRP
benchmark (200--480 customers). For related methodological context, see:

F. Li, B. Golden, and E. Wasil, "Very large-scale vehicle routing: new test
problems, algorithms, and results," *Computers & Operations Research*,
32(5), 1165--1179, 2005.
https://doi.org/10.1016/j.cor.2003.10.002
