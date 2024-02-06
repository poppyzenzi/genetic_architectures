require(GenomicSEM)
library(GenomicSEM)
library(readr)
# munged sum stats in "/exports/igmm/eddie/GenScotDepression/users/poppy/PRS/GWAS/"

#run multivariable LDSC to create the S and V matrices
psycho <- readRDS("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/ldsc/august_ipsych/LDSCoutput.rds")

#smooth the S matrix for EFA using the nearPD function in the Matrix package.
#require(Matrix)
#Ssmooth<-as.matrix((nearPD(psycho$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
#require(stats)
#EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")

#print the loadings
#EFA$loadings

#Specify the Genomic confirmatory factor model
model <- "F1 =~ NA * ANX + NEU + MDD
          F2 =~ NA * BIP + SCZ
          F3 =~ NA * ADHD + ASD

F1~~0*F2
F2~~0*F3
F3~~0*F1

P~~1*P

P =~ NA * F1 + F2 + F3

MDD ~~ c*MDD
c > 0.001"

#run the model
Psycho <- usermodel(psycho, estimation = "DWLS", model = model, 
            CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
Psycho
setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/efa_cfa/august_ipsych")
write.table(Psycho$results, file = "high_modelresult.txt", sep = "\t", row.names = FALSE)


