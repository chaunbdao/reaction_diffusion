# reactiondiffusion
MATLAB scripts for 1D stochastic reaction-diffusion simulations.

## Active Files

- `gen_pd_point.m`: serial calculation for one `(t_corr, t_v)` point; returns average valid patch spatial size.
- `gen_pd_point_par.m`: parallel version using `parfor` over independent runs; returns average valid patch spatial size, percent of runs with valid patches, and average patch lifetime.
- `rdsim.m`: phase-diagram driver over a `t_corr` by `t_v` grid; calls `gen_pd_point_par` and writes `pdmat`, `prpmat`, and `pltmat` CSV files.
- `periodic_test.m`: diagnostic script for monitoring the system behavior for a single run, test with different noise sources.
- `CC2periodic.m`: helper for merging connected components across periodic boundaries.
- `naive_rd_job.slurm.sh`: Rivanna Slurm script for running `rdsim`.

## Rivanna

Submit from the repo directory:

```bash
sbatch naive_rd_job.slurm.sh
```

The Slurm script requests one node, sets `numWorkers = SLURM_NTASKS - 1`, and passes that worker count into `rdsim`.

## Outputs

`rdsim` writes three parameter-labeled CSV files:

- `pdmat_...csv`: average valid patch spatial size.
- `prpmat_...csv`: percent of runs with at least one valid patch.
- `pltmat_...csv`: average valid patch lifetime in simulation time units.

## Deprecated

Older scripts are kept in `deprecated/` and are not part of the current phase-diagram workflow.

Generated CSV files are ignored by Git.
