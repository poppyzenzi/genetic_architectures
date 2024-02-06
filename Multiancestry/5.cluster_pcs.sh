#$ -N cluster_pcs
#$ -l h_vmem=3G
#$ -pe sharedmem 6
#$ -l h_rt=1:00:00
#$ -e logs
#$ -o logs
#$ -cwd

# Template for making PCs within an ancestry similarity group


# regions to exclude 
# - no MHC (6:25-35Mb)
# - no Chr.8 inversion (8:7-13Mb)
cat /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19.bim | awk '($1 == 6 && 25000000 <= $4 && $4 <= 35000000) || ($1 == 8 && 7000000 <= $4 && $4 <= 13000000) {print $2}' > data/regions.snplist

# list of strand-ambiguous SNPs to exclude
cat /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19.bim | awk '($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G") || ($5 == "A" && $6 == "T") ||  ($5 == "T" && $6 == "A") {print $2}' > data/strands.snplist

# prune SNPs for EAS sample
/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--bfile /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19 \
--exclude data/regions.snplist data/strands.snplist \
--snps-only 'just-acgt' \
--rm-dup 'exclude-all' \
--maf 0.05 \
--geno 0.02 \
--max-alleles 2 \
--hwe 1e-3 \
--indep-pairwise 200 100 0.2 \
--threads 6 \
--memory 17408 \
--keep-col-match-name 'ancestry_cluster' \
--keep-col-match data/ukb.randomforest.ancestries.tsv 'EAS' \
--out data/EAS

# Calculate PCs in EAS sample
/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--bfile /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19.bim \
--pca \
--threads 6 \
--memory 17408 \
--extract data/EAS.prune.in \
--keep-col-match-name 'ancestry_cluster' \
--keep-col-match data/ukb.randomforest.ancestries.tsv 'EAS' \
--out data/EAS_pcs
