#!/bin/bash

#SBATCH -p standard
#SBATCH --time=5-00:00:00
#SBATCH --account=cherngroup
#SBATCH --mail-type=end
#SBATCH --mail-user=chaunbdao@ucla.edu
#SBATCH --job-name=reactionDiffusionSimulations3
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
tv=20;                   # fixed tau_v
tcorr_start=0;           # start for correlation time sweep
tcorr_end=0.9;           # end for correlation time sweep (step = 0.1)
A_start=0.01;            # start for fixed-strength noise amplitude sweep
A_end=0.5;               # end for noise amplitude sweep (step = 0.01)

nRuns=50;                # number of runs collected for a given point
A_points=50;             # points for noise amplitude sweep
tcorr_points=10;         # points for correlation time sweep

# Run Matlab program
matlab -batch "rdsim3(${points},${physical_runtime},${nRuns},${tcorr_start},${tcorr_end},${tcorr_points},${A_start},${A_end},${A_points},${numWorkers},${dt},${tv}); exit;"
