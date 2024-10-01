#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 00_soupx
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 10:00:00
#SBATCH --output /projects/p31535/boles/img_scfrp_pilot/logs/%x.oe%j.log
#SBATCH --verbose

module load R/4.2.3

Rscript /projects/p31535/boles/img_scfrp_pilot/r_scripts/00_soupx.R