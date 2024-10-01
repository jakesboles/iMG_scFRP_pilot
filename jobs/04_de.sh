#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name 04_de
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 10:00:00
#SBATCH --output /projects/p31535/boles/img_scfrp_pilot/logs/%x.oe%j.log
#SBATCH --verbose

module load R/4.2.3

Rscript /projects/p31535/boles/img_scfrp_pilot/r_scripts/04_de.R