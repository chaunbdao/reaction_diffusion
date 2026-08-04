# reactiondiffusion
MATLAB scripts for 1D stochastic reaction-diffusion simulations.

## Active Files

- `gen_pd_point.m`: serial calculation for one `(t_corr, t_v)` point; returns average valid patch spatial size.
- `gen_pd_point_par.m`: parallel version using `parfor` over independent runs; returns average valid patch spatial size, percent of runs with valid patches, and average patch lifetime.
- `gen_pd_point_par2.m`: optimized worker used by `rdsim2`, `rdsim3`, and `rdsim4`; precomputes timestep matrices, reuses a sparse LU solver per worker, and analyzes snapshots every `0.1` time units. Optional arguments select the noise mode and initial inhibitor value.
- `rdsim.m`: phase-diagram driver over a `t_corr` by `t_v` grid; calls `gen_pd_point_par` and writes `pdmat`, `prpmat`, and `pltmat` CSV files.
- `rdsim2.m`: physical-runtime phase-diagram driver with supplied `dt` and noise strength `A`.
- `rdsim3.m`: physical-runtime phase-diagram driver over fixed-strength noise amplitude `A` and `t_corr`, with a fixed `t_v`.
- `rdsim4.m`: phase-diagram driver over `t_corr` and `t_v` using fixed instantaneous variance, `Y(0) = -6`, and a supplied fixed `sigma_eta`.
- `periodic_test_optimized.m`: main single-run diagnostic; uses precomputed timestep matrices and a reused sparse LU decomposition.
- `CC2periodic.m`: helper for merging connected components across periodic boundaries.
- `naive_rd_job.slurm.sh`: Rivanna Slurm script for running `rdsim`.
- `naive_rd_job2.slurm.sh`: Rivanna Slurm script for running `rdsim2`.
- `naive_rd_job3.slurm.sh`: Rivanna Slurm script for running `rdsim3`.
- `rdslurm_job4.slurm.sh`: Rivanna Slurm script for running `rdsim4` with `sigma_eta = 0.5`.

## Rivanna

Submit from the repo directory:

```bash
sbatch naive_rd_job.slurm.sh
sbatch naive_rd_job2.slurm.sh
sbatch naive_rd_job3.slurm.sh
sbatch rdslurm_job4.slurm.sh
```

`naive_rd_job.slurm.sh` sets `numWorkers = SLURM_NTASKS - 1` for `rdsim`.

`naive_rd_job2.slurm.sh` requests one MATLAB task with six CPUs and sets `numWorkers = SLURM_CPUS_PER_TASK - 1` for `rdsim2`.

`naive_rd_job3.slurm.sh` uses the same allocation pattern for `rdsim3`.

`rdsim2` keeps the original sweep arguments and adds `dt` and `A`:

```matlab
rdsim2(points,physical_runtime,nruns,tcorr_start,tcorr_end,tcorr_points,tv_start,tv_end,tv_points,workers,dt,A)
```

`rdsim3` sweeps `A` and `t_corr` while holding `t_v` fixed:

```matlab
rdsim3(points,physical_runtime,nruns,tcorr_start,tcorr_end,tcorr_points,A_start,A_end,A_points,workers,dt,tv)
```

`rdsim4` sweeps `t_corr` and `t_v` with fixed-instantaneous-variance noise:

```matlab
rdsim4(points,physical_runtime,nruns,tcorr_start,tcorr_end,tcorr_points,tv_start,tv_end,tv_points,workers,dt,sigma_eta)
```

## Outputs

`rdsim` writes three parameter-labeled CSV files:

- `pdmat_...csv`: average valid patch spatial size.
- `prpmat_...csv`: percent of runs with at least one valid patch.
- `pltmat_...csv`: average valid patch lifetime in simulation time units.

`rdsim2` writes `pdmat2`, `prpmat2`, and `pltmat2` files labeled with `dt`, physical runtime, `A`, and snapshot interval.

`rdsim3` writes `pdmat3`, `prpmat3`, and `pltmat3` files. Rows correspond to `A`, columns correspond to `t_corr`, and filenames include the fixed `t_v`. When `t_corr = 0`, `gen_pd_point_par2` uses the corresponding Gaussian white-noise limit.

`rdsim4` writes `pdmat4`, `prpmat4`, and `pltmat4` files. Rows correspond to `t_v`, columns correspond to `t_corr`, and filenames include `active_fixed_variance`, `sigma_eta`, and the initial `Y` value.

## Collected Data

`collected_data/` contains phase-diagram CSV blocks for `t_v = 1:50` and `t_corr = 0.5:0.1:5.4`.

- `plottertest.ipynb`: combines the CSV blocks and plots phase diagrams and selected `t_corr` slices.

## Deprecated

Older scripts are kept in `deprecated/` and are not part of the current phase-diagram workflow.

This includes the original `periodic_test.m` diagnostic.

Generated CSV files are ignored by Git.
