require(GenomicSEM)
library(GenomicSEM)
library(readr)

# read in covstruct and sumstats
# using mdd ispsych - check files are from these subdirs
cat("Reading covstruct and sumstats...\n")
covstruct <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/august_ipsych/LDSCoutput.rds")
sumstats <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/common_factor/august_ipsych/PSYCH_sumstats.rds")

# specify model with SNP effects
cat("Specifying model...\n")
model <- "F1 =~ ANX + NEU + MDD 
	  F2 =~ a*BIP + a*SCZ 
	  F3 =~ b*ADHD + b*ASD 
F1~~F2
F2~~F3
F1~~F3

F1~SNP
F2~SNP
F3~SNP

MDD~~c*MDD
c>0.001"

# set wd for saving tsvs
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/corr_factor/")

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST"))

cat(paste("Running task", k, " of ", max_tasks, "\n"))

out_tsv <- paste0('gwas/mood/F1_', k, '.tsv')
out_tsv2 <- paste0('gwas/psychotic/F2_', k, '.tsv')
out_tsv3 <- paste0('gwas/neurodev/F3_', k, '.tsv')

dir.create(dirname(out_tsv), showWarnings=FALSE)
dir.create(dirname(out_tsv2), showWarnings=FALSE)
dir.create(dirname(out_tsv3), showWarnings=FALSE)

# get subset of SNPs to analyze in this task
nsnps <- nrow(sumstats)
snps_per_task <- ceiling(nsnps / max_tasks)


# index of task IDs repeated for number of SNPs per ask, clip the last task# to the actual number of snps
task_allocations <- rep(1:max_tasks, each=snps_per_task)[seq.int(nsnps)]
SNPs <- sumstats[which(task_allocations == k),]
rm(sumstats)


# get SNP effects on each mood parcel
cat("Running GWAS...\n")
gwas <- userGWAS(covstruc=covstruct, toler = 1e-50, model = model, 
                                        SNPs=SNPs,
					sub=c("F1~SNP", "F2~SNP","F3~SNP"),
                                        smooth_check=TRUE,
                                        parallel=TRUE, cores=8)

# Save the complete output of userGWAS to a file
cat("Saving the complete output of GWAS...\n")
saveRDS(gwas, "complete_gwas_output.rds")

# Write specific components to TSV files
cat("Writing results to files...\n")
write_tsv(gwas[[1]], out_tsv)
write_tsv(gwas[[2]], out_tsv2)
write_tsv(gwas[[3]], out_tsv3)


