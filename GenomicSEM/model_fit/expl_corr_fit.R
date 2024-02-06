require(GenomicSEM)
library(GenomicSEM)
library(readr)
# munged sum stats in "/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/"
#run multivariable LDSC to create the S and V matrices
psycho <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/august_ipsych/LDSCoutput.rds")

#smooth the S matrix for EFA using the nearPD function in the Matrix package.
require(Matrix)
Ssmooth<-as.matrix((nearPD(psycho$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")
#print the loadings
EFA$loadings

# confirmatory
# correlated factor model
#Specify the Genomic confirmatory factor model

model <- "F1 =~ ANX + NEU + MDD 
          F2 =~ a*BIP + a*SCZ
          F3 =~ b*ADHD + b*ASD
F1~~F2
F2~~F3
F1~~F3


MDD~~c*MDD
c>0.001"


#run the model
Psycho <- usermodel(psycho, estimation = "DWLS", model = model, 
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Psycho

setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/efa_cfa/august_ipsych")
write.table(Psycho$results, file = "corr_modelresult.txt", sep = "\t", row.names = FALSE)


