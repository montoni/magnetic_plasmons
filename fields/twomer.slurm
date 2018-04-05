#!/bin/bash 
## Job Name 
#SBATCH --job-name=scafield
## Allocation Definition
#SBATCH --account=chem
#SBATCH --partition=chem
## Resources 
## Nodes 
#SBATCH --nodes=1       
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=28 
## Walltime (ten minutes) 
#SBATCH --time=48:00:00 
## Memory per node 
#SBATCH --mem=16G 
## Specify the working directory for this job 
#SBATCH --workdir=/gscratch/chem/montoni/magnets_narrow/fields 

module load matlab_2017a

matlab -nodisplay -r twomer_NS
