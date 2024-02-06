require(GenomicSEM)
library(GenomicSEM)
library(readr)

# read in covstruct and sumstats

LDSCoutput <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/august_ipsych/LDSCoutput.rds")
PSYCH_sumstats <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/common_factor/august_ipsych/PSYCH_sumstats.rds")


covstruct <- LDSCoutput
sumstats <- PSYCH_sumstats

# specify model with SNP effects

model <- "F1 =~ ANX + NEU + MDD 
	  F2 =~ BIP + SCZ 
	  F3 =~ ADHD + ASD 

P =~ F1 + F2 + F3

P~SNP

MDD~~c*MDD
c>0.001"


# set wd for saving tsvs
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/high_factor/august_ipsych")

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST"))

cat(paste("Running task", k, " of ", max_tasks, "\n"))
out_tsv <- paste0('gwas/', k, '.tsv')
dir.create(dirname(out_tsv), showWarnings=FALSE)


# get subset of SNPs to analyze in this task
nsnps <- nrow(sumstats)
snps_per_task <- ceiling(nsnps / max_tasks)


# index of task IDs repeated for number of SNPs per ask, clip the last task# to the actual number of snps
task_allocations <- rep(1:max_tasks, each=snps_per_task)[seq.int(nsnps)]
SNPs <- sumstats[which(task_allocations == k),]
rm(sumstats)


gwas <- userGWAS(covstruc=covstruct, toler = 1e-40, model = model, 
                                        SNPs=SNPs,
					sub=c("P~SNP"),
                                        smooth_check=TRUE,
                                        parallel=TRUE, cores=8)
# save first item (P~SNP)
write_tsv(gwas[[1]], out_tsv)
