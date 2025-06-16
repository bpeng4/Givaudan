#!/bin/bash
#SBATCH --time=168:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=5G       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=Trim
#SBATCH --error=/work/benson/bpeng4/Givaudan/err/job.%J.err
#SBATCH --output=/work/benson/bpeng4/Givaudan/err/job.%J.output
#SBATCH --partition=benson,batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36



# Load the anoconda environment
module purge
ml anaconda
conda activate mydada2

Rscript Trim.R
