#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ben.sherwood@ku.edu
#SBATCH --partition=sherwood
#SBATCH --time=21-00:00:00
#SBATCH --mem=10000
module load R 
R CMD BATCH --vanilla AsymWild.r