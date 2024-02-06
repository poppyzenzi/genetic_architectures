require(GenomicSEM)
library(GenomicSEM)
library(readr)
# munged sum stats in "/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/"


#run multivariable LDSC to create the S and V matrices
psycho <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/august_ipsych/LDSCoutput.rds")


#Specify the Genomic confirmatory factor model
model <- 'F1 =~ ANX + NEU + MDD + BIP + SCZ + ADHD + ASD

MDD ~~ c*MDD
c > 0.001'

#run the model
Psycho <- usermodel(psycho, estimation = "DWLS", model = model, 
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

#print the results
Psycho

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/efa_cfa/august_ipsych/")

write.table(Psycho$results, file = "cf_modelresult.txt", sep = "\t", row.names = FALSE)


