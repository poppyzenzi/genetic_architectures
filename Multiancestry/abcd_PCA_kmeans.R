#ABCD3.0 - SELECTION OF WHITE ONLY POP USING KMEANS CLUSTERING OF WHOLESAMPLE GENETIC PCs
library(tidyverse)
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(patchwork)

#import ethnicity reported 
demo <- read.csv("abcd_p_demo.csv", header=TRUE) 
demo <- demo[!duplicated(demo$src_subject_id), ]

#import PCs generated from whole sample PCA
setwd("/exports/eddie/scratch/s2421111/abcd/genetics_QCed/genotyped/")
pca <- read.delim("ABCD3.0_wholesample_genotyped_PCA.eigenvec", header=F, sep="")
colnames(pca)[1] <- "FID"
colnames(pca)[2] <- "src_subject_id"
colnames(pca)[3:17] <- rep(1:15)
colnames(pca)[3:17] <-paste0('PC_',colnames(pca)[3:17])

#fix the 2IDs that are incorrect  -> `NDAR_INVF3FYXH1G and NDARINVPWLFYWTX
pca$src_subject_id[grep("`NDAR_INVF3FYXH1G",pca$src_subject_id)] <- "NDAR_INVF3FYXH1G"
pca$src_subject_id[grep("NDARINVPWLFYWTX",pca$src_subject_id)] <- "NDAR_INVPWLFYWTX"

#check scree plot
eigenval <- read.delim("ABCD3.0_wholesample_genotyped_PCA.eigenval", header=F, sep="") %>%  rownames_to_column()
colnames(eigenval) <- c("PC","eigenval")

eigenval %>%
  ggplot(aes(x=as.numeric(fct_inorder(PC)),y=eigenval))+
  geom_col()+
  geom_line()+
  theme_minimal()

#K-MEANS CLUSTERING of PCs obtained from plink whole sample ====
#follow steps: https://static-content.springer.com/esm/art%3A10.1038%2Fng.3768/MediaObjects/41588_2017_BFng3768_MOESM40_ESM.pdf
set.seed(2021)
# k means for pca cols 8 pc cols starting from i+2
kmeans_results <- list()
for (i in 1:8) {
        km_result <- kmeans(pca[, i+2], 4, nstart=10)         
	kmeans_results[[paste0("km.pc",i)]] <- km_result }

 
# df of id cols, pca cols, cluster assignments
pca.cluster <- data.frame(pca[, 1:2], pca[, 3:10], stringsAsFactors = FALSE)
for (i in 1:8) {
   cluster_assignments <- kmeans_results[[paste0("km.pc", i)]]$cluster
   pca.cluster[[paste0("km.pc", i)]] <- cluster_assignments
 }
 
# get self/parent(?) reported ethnicities
race.reported <- data.frame(src_subject_id= demo$src_subject_id,
                             race = demo$race_ethnicity) #white, black, hispanic, asians,others in order
# merge report with df
pca.race.cluster <- merge(pca.cluster,race.reported,by="src_subject_id",all.x=T)
head(pca.race.cluster)

# freq tables  with reported race
for (i in 1:8) {
   freq_table <- table(pca.race.cluster$race, pca.race.cluster[[paste0("km.pc", i)]])
   print(freq_table)
 } 


setwd("/exports/igmm/eddie/GenScotDepression/users/poppy/abcd/multian/")

options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 300)
subplot_list <- list()

# Create four subplots of 2 PCs each
for (i in 1:4) {
  p <- ggplot(pca.race.cluster, aes(x = .data[[paste0("PC_", (i - 1) * 2 + 1)]],
                                    y = .data[[paste0("PC_", (i - 1) * 2 + 2)]],
                                    color = race)) +
    geom_point(size=0.5) +
    labs(x = paste0("PC_", (i - 1) * 2 + 1), y = paste0("PC_", (i - 1) * 2 + 2)) +
    theme_minimal() +
    theme(legend.position = "none", aspect.ratio = 1)
  subplot_list[[i]] <- p
}

# Combine the subplots into a grid and add a common legend
subplot_grid <- wrap_plots(plotlist = subplot_list, ncol = 2) +
  plot_layout(guides = "collect")
 
# Save the plot to a file with improved resolution
ggsave("pca_subplots.png", subplot_grid, dpi = 300)



#write.table(black.only, "abcd4.0_blackonly_genotyped.txt", sep = ' ', quote=F, row.names=F, col.names=F)

pca.cluster$cluster.label <- ifelse(pca.cluster$km.pc1 == 3 & pca.cluster$km.pc2 == 2, 'european', 
			ifelse(pca.cluster$km.pc1 ==2 & pca.cluster$km.pc2 ==2, 'african', 
			ifelse(pca.cluster$km.pc1 ==1 & pca.cluster$km.pc2 ==4, 'asian','other')))

p <- ggplot(pca.cluster, aes(x=PC_1, y=PC_2, color=cluster.label)) + 
	geom_point() + labs(x = "PC1", y = "PC2") + ggtitle("PCA Plot") + theme_minimal()
#png("pca_plot2.png")
#print(p)
#dev.off()

#african.only <- pca.cluster %>% filter(cluster.label == 'african') %>%  select(FID, src_subject_id)
#european.only <- pca.cluster %>% filter(cluster.label == 'european') %>% select(FID, src_subject_id)



