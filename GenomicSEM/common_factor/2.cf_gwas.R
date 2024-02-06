# Common factor GWAS

require(GenomicSEM)
library(GenomicSEM) 
library(readr) 

# read in covstruct and sumstats
setwd('/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/common_factor/LOO/')
LDSCoutput <- readRDS("LDSCoutput_no_asd.rds")
PSYCH_sumstats <- readRDS("LOO_asd_sumstats.rds")

covstruct <- LDSCoutput 
sumstats <- PSYCH_sumstats 

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST")) 

cat(paste("Running task", k, " of ", max_tasks, "\n")) 

out_tsv <- paste0('gwas/no_asd/', k, '.tsv') 

dir.create(dirname(out_tsv), showWarnings=FALSE) 

# get subset of SNPs to analyze in this task
nsnps <- nrow(sumstats) 
snps_per_task <- ceiling(nsnps / max_tasks) 

# index of task IDs repeated for number of SNPs per ask, clip the last task# to the actual number of snps
task_allocations <- rep(1:max_tasks, each=snps_per_task)[seq.int(nsnps)] 
SNPs <- sumstats[which(task_allocations == k),] 
rm(sumstats) 

gwas <- commonfactorGWAS(covstruc=covstruct, 
					SNPs=SNPs, 
					smooth_check=TRUE, 
					parallel=TRUE, cores=12) 
write_tsv(gwas, out_tsv)
