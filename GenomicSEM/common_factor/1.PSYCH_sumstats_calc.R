# step 3: prep sum stats for GWAS
# model with SNP effects = multivariate GWAS sumstats

#require(GenomicSEM)
library(GenomicSEM)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/")

files<-c(
	"adhd_23.gz",
	"pgc_bip_qcd.txt",
	#"pgc_scz_qcd.txt",
        "iPSYCH-PGC_ASD_Nov2017_INFO0.8_nodup_noambig.gz",
	"mdd2023ipsych_with_N.gz",
	"luc_neur_ss.txt",
	"ANX_withNeff_nodup_noambig.gz")

ref= "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/reference.1000G.maf.0.005.txt"

trait.names<-c(
	'ADHD', 
	'BIP', 
	'SCZ', 
	'ASD',
	'MDD',
	'NEU' 
	#'ANX'
)

# whether SE is logistic - check README file
se.logit=c(
	#T,
	T,T,T,T,T,T)

# Whether the phenotype was a dichotomous outcome for which there are *only Z-statistics* in the
# summary statistics file -or- it was a dichotomous outcome analyzed using an OLS estimator as in UKB for some using hail software
linprob=NULL
# Whether the phenotype was a continuous outcome analyzed using an observed least square (OLS; i.e., linear)
OLS=NULL

info.filter=0.8
maf.filter=0.01

# if ss file has sample size col then can put NA, make sure NEff for dichiotomous traits

N=c(
	103135.527,
	NA,
	NA,
	43777.81914,
	NA,
	329821
	#NA
)

# cont trait, if MAF col not available
betas=NULL

# sum stats -> use to generate p-factor PRS
PSYCH_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,
                            se.logit=se.logit,OLS=OLS,linprob=linprob,N=N,
                            betas=betas,info.filter=info.filter,maf.filter=maf.filter,
                            keep.indel=FALSE,parallel=FALSE,cores=NULL)

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/common_factor/LOO")

save(PSYCH_sumstats, file = "LOO_anx_sumstats.RData")
saveRDS(PSYCH_sumstats, file = "LOO_anx_sumstats.rds")

# =====================================================================================================================

