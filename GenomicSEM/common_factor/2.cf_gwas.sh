#!/bin/bash
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 12
#$ -cwd
#$ -N comfac_no_asd
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m as
#$ -t 1-100

# load modules

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0 


# if running this script put task numbers in shell commands
Rscript 2.cf_gwas.R
