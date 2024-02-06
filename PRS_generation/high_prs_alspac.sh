#!/bin/bash
#$ -N alspac_high_factor_prs
#$ -l h_rt=1:00:00
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/gsem/job_logs
#$ -M p.grimes@ed.ac.uk
#$ -m baes

. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0

# base data and saving output
HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/high_factor
# target data
ALSPAC=/exports/igmm/eddie/GenScotDepression/users/poppy/alspac/prs

# alspac, check time point
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $HOME/gwas/gwastable.txt \
        --target $ALSPAC/data_QC_filtered \
        --beta \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat est --pvalue Pval_Estimate \
        --pheno-file $ALSPAC/pheno_file_alspac_all.txt \
        --pheno-col t1,t2,t3,t4 \
        --binary-target F,F,F,F \
	--perm 10000 \
        --fastscore \
	--ignore-fid \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
        --out $HOME/prs/alspac/alspac_highfac_prs$(date +%m%d)


