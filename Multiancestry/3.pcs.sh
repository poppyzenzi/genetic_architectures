#$ -N pcs
#$ -hold_jid tkg
#$ -l h_vmem=3G
#$ -pe sharedmem 6
#$ -l h_rt=1:00:00
#$ -e logs
#$ -o logs
#$ -cwd

# extract pruned SNPs from ABCD
/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--bfile /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19 \
--extract data/1kg_ld.prune.in \
--make-bed \
--out data/background \
--threads 6 \
--memory 17408 

# extract pruned SNPs from 1KG file
/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile /exports/igmm/eddie/GenScotDepression/data/resources/1000g/PGEN/all_phase3 'vzs' \
--extract data/1kg_ld.prune.in \
--make-bed \
--max-alleles 2 \
--out data/1kg \
--threads 6 \
--memory 17408

touch data/background_1kg_check.missnp

# check merge with 1KG
/exports/igmm/eddie/GenScotDepression/local/bin/plink \
--bfile data/1kg \
--bmerge data/background \
--merge-mode 6 \
--out data/background_1kg_check \
--threads 6 \
--memory 17408

# remove mismatched SNPs
/exports/igmm/eddie/GenScotDepression/local/bin/plink \
--bfile data/1kg \
--exclude data/background_1kg_check.missnp \
--make-bed \
--out data/1kg_checked \
--threads 6 \
--memory 17408

/exports/igmm/eddie/GenScotDepression/local/bin/plink \
--bfile data/background \
--exclude data/background_1kg_check.missnp \
--make-bed \
--out data/background_checked \
--threads 6 \
--memory 17408

# merge checked files
/exports/igmm/eddie/GenScotDepression/local/bin/plink \
--bfile data/1kg_checked \
--bmerge data/background_checked \
--out data/background_1kg \
--threads 6 \
--memory 17408

# project combined sample into 1KG PCs
/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--bfile data/background_1kg \
--read-freq data/1kg_pca.acount \
--score data/1kg_pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-15 \
--out data/background_1kg_projection \
--threads 6 \
--memory 17408


