#$ -N tkg
#$ -l h_vmem=3G
#$ -pe sharedmem 6
#$ -l h_rt=1:00:00
#$ -e logs
#$ -o logs
#$ -cwd

## Calculate PCS for 1000 Genomes
# 1000 is the reference data file 
mkdir -p data

# regions to exclude (from abcd imputed, but check this?)
# - no MHC (6:25-35Mb)
# - no Chr.8 inversion (8:7-13Mb)
cat /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19.bim | awk '($1 == 6 && 25000000 <= $4 && $4 <= 35000000) || ($1 == 8 && 7000000 <= $4 && $4 <= 13000000) {print $2}' > data/regions.snplist


# list of strand-ambiguous SNPs to exclude
cat /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19.bim | awk '($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G") || ($5 == "A" && $6 == "T") ||  ($5 == "T" && $6 == "A") {print $2}' > data/strands.snplist


# calculate PCs in 1KG sample
/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile /exports/igmm/eddie/GenScotDepression/data/resources/1000g/PGEN/all_phase3 'vzs' \
--extract /exports/eddie/scratch/s2421111/genetics_QCed/imputed/hg19/wholesample/ABCD3.0_imputed_wholesample_hg19.bim \
--exclude data/regions.snplist data/strands.snplist \
--snps-only 'just-acgt' \
--rm-dup 'exclude-all' \
--maf 0.05 \
--geno 0.02 \
--max-alleles 2 \
--hwe 1e-3 \
--indep-pairwise 200 100 0.2 \
--out data/1kg_ld \
--threads 6 \
--memory 17408 

/exports/igmm/eddie/GenScotDepression/local/bin/plink2 \
--pfile /exports/igmm/eddie/GenScotDepression/data/resources/1000g/PGEN/all_phase3 'vzs' \
--extract data/1kg_ld.prune.in \
--freq counts \
--pca allele-wts \
--out data/1kg_pca \
--threads 6 \
--memory 17408

