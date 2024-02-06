require(GenomicSEM)
library(GenomicSEM)
library(readr)
#run multivariable LDSC to create the S and V matrices
psycho <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/august_ipsych/LDSCoutput.rds")

#Specify the Genomic confirmatory factor model
# bifactor model 
model <- "F1 =~ NA*ANX + NEU + MDD
          F2 =~ NA*BIP + SCZ
          F3 =~ NA*ADHD + ASD

P =~ NA*ANX + NEU + MDD + BIP + SCZ + ADHD + ASD

F1 ~ SNP
F2 ~ SNP
F3 ~ SNP
P ~ SNP

F1~~1*F1
F2~~1*F2
F3~~1*F3

P~~0*F1
P~~0*F2
P~~0*F3

F1~~0*F2
F2~~0*F3
F3~~0*F1

MDD ~~ c*MDD
c > 0.001"


#run the model
Psycho <- usermodel(psycho, estimation = "DWLS", model = model, sub=c("P~SNP"),
 CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Psycho
Psycho$modelfit
Psycho$results

