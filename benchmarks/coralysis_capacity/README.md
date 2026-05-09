# Coralysis Capacity Benchmark

This benchmark measures `sn_run_cluster()` runtime and peak resident memory for
the Coralysis and Coralysis2 integration backends on downsampled or resampled
PBMC Seurat objects.

The driver creates benchmark input objects from
`/home/sduan/projects/immune-atlas/data/processed/pbmc.qs`, relabels the
`sample` column into a requested number of synthetic batches, and runs each
configuration in a fresh R process under `/usr/bin/time -v`.

Example quick pilot:

```sh
Rscript benchmarks/coralysis_capacity/run_grid.R \
  --cell_counts 3000,10000 \
  --batch_counts 3,10 \
  --methods coralysis,coralysis2 \
  --threads 1 \
  --L 2 \
  --k 4 \
  --train_with_bnn false \
  --outdir benchmarks/coralysis_capacity/results_quick
```

Example larger Coralysis2-only run:

```sh
Rscript benchmarks/coralysis_capacity/run_grid.R \
  --cell_counts 30000,60000,120000 \
  --batch_counts 10,30 \
  --methods coralysis2 \
  --threads 1,4 \
  --L 3 \
  --k 8 \
  --outdir benchmarks/coralysis_capacity/results_large
```

Outputs:

- `summary.csv`: one row per run, including exit status, elapsed time, peak RSS,
  object dimensions, method, thread count, and error text if the run failed.
- `inputs/`: generated benchmark Seurat objects.
- `runs/*.stdout.txt`, `runs/*.stderr.txt`, `runs/*.time.txt`, `runs/*.json`:
  per-run logs and metrics.

The synthetic batch relabeling isolates scaling with batch count. It is not a
biological benchmark of sample heterogeneity.
