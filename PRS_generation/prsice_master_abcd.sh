#!/bin/bash
#$ -N prsice_abcd_all_prs
#$ -l h_rt=12:00:00
#$ -l h_vmem=128G
#$ -pe sharedmem 1
#$ -e /exports/igmm/eddie/GenScotDepression/users/poppy/abcd/job_logs
#$ -o /exports/igmm/eddie/GenScotDepression/users/poppy/abcd/job_logs
#$ -M s2421111@ed.ac.uk
#$ -m baes


# load modules
. /etc/profile.d/modules.sh
module load igmm/apps/PRSice/2.1.11
module unload igmm/apps/R/3.5.1
module load igmm/apps/R/4.1.0
 

HOME=/exports/igmm/eddie/GenScotDepression/users/poppy/abcd
BASE=/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS


sumstats_files=(
    "$BASE/mdd3_ss_MAF.txt"
    "$BASE/adhd_23.gz"
    "$BASE/iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz"
    "$BASE/ANX_withNeff_nodup_noambig.gz"
    "$BASE/pgc_scz_qcd.txt"
    "$BASE/pgc_bip_qcd.txt"
    "$BASE/luciano_neur_ss.txt"
    "$BASE/mdd3_ss_MAF.txt"
)

###### before running need to make sure columns in base files have same headers (ID, CHR, POS, A1 etc...)
for sumstats in "${sumstats_files[@]}"; do
    Rscript /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice.R \
        --prsice /exports/igmm/software/pkg/el7/apps/PRSice/2.1.11/PRSice_linux \
        --base "$sumstats" \
        --target "$HOME/ABCD3.0_imputed_whiteonly_MAF0.01_unrelated" \
        --binary-target F,F,F,F,F,F,F,F \
        --pheno-file "$HOME/phenotype/pheno_file_bpm_all.txt" \
        --snp ID --chr CHR --bp POS --A1 A1 --A2 A2 --stat BETA --pvalue PVAL \
        --beta \
        --perm 10000 \
        --pheno-col t0.5,t1,t1.5,t2,t2.5,t3,3.5,t4 \
        --ignore-fid \
        --all-score \
        --print-snp \
        --fastscore \
        --bar-levels 0.00000005,0.000001,0.0001,0.001,0.01,0.05,0.1,0.2,0.5,1 \
	--out "$HOME/output_aug_loop/${filename_without_extension}_$(date +%m%d)"
done

