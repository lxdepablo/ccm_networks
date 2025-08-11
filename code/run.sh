#!/bin/bash
#SBATCH -p atesting # Partition or queue
#SBATCH --qos=testing
#SBATCH --job-name=generate_data # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luis.depablo@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=8 # Run on n CPUs
#SBATCH --mem=16gb # Memory limit
#SBATCH --time=0:30:00 # Time limit hrs:min:sec
#SBATCH --output=log.out # Standard output and error log
#SBATCH --error=log.err # %j inserts job number

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
date

module purge
module load miniforge
mamba activate ccm_env

Rscript generate_xmaps.R
Rscript xmap_analysis.R
