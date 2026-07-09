# reactiondiffusion
MATLAB scripts for 1D stochastic reaction-diffusion simulations.

## Active Files

- `gen_pd_point.m`: serial calculation for one `(t_corr, t_v)` point; returns average valid patch spatial size.
- `gen_pd_point_par.m`: parallel version using `parfor` over independent runs; requires an active MATLAB parallel pool.
- `rdsim.m`: phase-diagram driver over a `t_corr` by `t_v` grid; calls `gen_pd_point_par` and writes `pdmat.csv`.
- `periodic_test.m`: diagnostic script for monitoring the system behavior for a single run, test with different noise sources.
- `CC2periodic.m`: helper for merging connected components across periodic boundaries.
- `naive_rd_job.slurm.sh`: Rivanna Slurm script for running `rdsim`.

## Rivanna

Submit from the repo directory:

```bash
sbatch naive_rd_job.slurm.sh
```

The Slurm script requests one node, sets `numWorkers = SLURM_NTASKS - 1`, and passes that worker count into `rdsim`.

## Deprecated

Older scripts are kept in `deprecated/` and are not part of the current phase-diagram workflow.

Generated CSV files are ignored by Git.
