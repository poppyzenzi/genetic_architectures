#!/bin/bash
#$ -N LDSC_no_asd_sumstats
#$ -l h_rt=1:00:00
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes
#$ -cwd

# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript LDSC_CFM.R
