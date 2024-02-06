#!/bin/bash
#$ -N high_factor_prs
#$ -l h_rt=2:00:00
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
ABCD=/exports/igmm/eddie/GenScotDepression/users/poppy/abcd

# abcd bpm
Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base $HOME/gwas/gwastable.txt \
        --target $ABCD/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated \
        --beta \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat est --pvalue Pval_Estimate \
        --pheno-file $ABCD/phenotype/pheno_file_bpm_all.txt \
        --pheno-col t0.5,t1,t1.5,t2,t2.5,t3,t3.5,t4 \
        --binary-target F,F,F,F,F,F,F,F \
	--perm 10000 \
	--ignore-fid \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
        --all-score  \
        --print-snp \
      --out $HOME/prs/abcd/abcd_highfac_prs$(date +%m%d)
