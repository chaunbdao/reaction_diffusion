#!/bin/bash

#SBATCH -p standard
#SBATCH --time=1-00:00:00
#SBATCH --account=cherngroup
#SBATCH --mail-type=end
#SBATCH --mail-user=chaunbdao@ucla.edu
#SBATCH --job-name=reactionDiffusionSimulations
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --chdir="/home/cnd4xt/chern_group/reaction_diffusion/reaction_diffusion"

# Load Matlab environment
module purge
module load matlab

export slurm_ID="${SLURM_JOB_ID}"

# Set workers to one less that number of tasks (leave 1 for master process)
export numWorkers=$((SLURM_NTASKS-1))


# Input parameters for Matlab function
points=1000; 		# number of points in the simulation
runtime=5000; 		# number of timesteps
tcorr_start=0.5;	# start for correlation time sweep
tcorr_end=5.5;		# end for correlation time sweep (inclusive)
tv_start=1;			# start for tau_v sweep
tv_end=50;			# end for tau_v sweep

nRuns=50; 			# number of collected for a given
tv_points=51;		# points for tau_v sweep
tcorr_points=51;	# points for correlation time sweep

# Run Matlab single core program
matlab -batch "rdsim(${points},${runtime},${nRuns},${tcorr_start},${tcorr_end},${tcorr_points},${tv_start},${tv_end},${tv_points},${numWorkers}); exit;"
