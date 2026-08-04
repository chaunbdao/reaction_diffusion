#!/bin/bash

#SBATCH -p standard
#SBATCH --time=5-00:00:00
#SBATCH --account=cherngroup
#SBATCH --mail-type=end
#SBATCH --mail-user=chaunbdao@ucla.edu
#SBATCH --job-name=reactionDiffusionFixedVariance4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --chdir="/home/cnd4xt/chern_group/reaction_diffusion/reaction_diffusion"

# Load Matlab environment
module purge
module load matlab

export slurm_ID="${SLURM_JOB_ID}"

# Set workers to one less than the allocated CPUs (leave 1 for the client).
export numWorkers=$((SLURM_CPUS_PER_TASK-1))


# Input parameters for Matlab function
points=1000;             # number of points in the simulation
physical_runtime=500;    # total physical simulation time
dt=0.01;                 # integration timestep
sigma_eta=0.5;           # fixed instantaneous noise standard deviation
tcorr_start=0.5;         # start for positive correlation-time sweep
tcorr_end=5.4;           # end for correlation-time sweep (inclusive)
tv_start=1;              # start for tau_v sweep
tv_end=50;               # end for tau_v sweep

nRuns=50;                # number of runs collected for a given point
tv_points=50;            # points for tau_v sweep
tcorr_points=50;         # points for correlation-time sweep

# Run Matlab program. rdsim4 fixes noise_mode to active_fixed_variance and Y(0) to -6.
matlab -batch "rdsim4(${points},${physical_runtime},${nRuns},${tcorr_start},${tcorr_end},${tcorr_points},${tv_start},${tv_end},${tv_points},${numWorkers},${dt},${sigma_eta}); exit;"
