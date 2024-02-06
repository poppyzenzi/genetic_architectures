# cluster ABCD participants with 1KG ancestries

library(readr)
library(dplyr)
library(ggplot2)
library(randomForest)

# projected PCs for combined sample
pcs <- read_tsv('data/background_1kg_projection.sscore')

# ethnic background in UKB
abcd_background <- read_tsv('data/ethnic_background.txt')

# 1KG population groups
pops <- read_tsv('/exports/igmm/eddie/GenScotDepression/data/resources/1000g/PGEN/all_phase3.psam')

superpops <- pops %>%
select(IID=`#IID`, ancestry=SuperPop)

# merge ethnic background of ABCD individuals
pcs_superpops <-
abcd_background %>%
filter(IID %in% pcs$IID) %>%
transmute(IID=as.character(IID), ancestry=ethnicity) %>%
bind_rows(superpops) %>%
left_join(pcs, by='IID')

dir.create('figures')

ggplot(pcs_superpops, aes(x=PC1_AVG, PC2_AVG, colour=ancestry)) +
geom_point()
ggsave('figures/superpops.pc12.png')

# training set of 1KG populations
training_pcs <- pcs_superpops %>%
    filter(ancestry %in% c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'))
    
rf <- randomForest(as.factor(ancestry) ~ PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG, data=training_pcs)

testing_pcs <- pcs_superpops %>%
    filter(!ancestry %in% c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'))
    
testing_pred <- predict(rf, testing_pcs)

testing_prob <- predict(rf, testing_pcs, type='prob')

ancestry_pred <- testing_pcs %>%
                 mutate(ancestry_hat=testing_pred,
                        ancestry_prob=apply(testing_prob, 1, max))
                        
ancestry_clusters <- ancestry_pred %>%
    mutate(ancestry_cluster=if_else(ancestry_prob >= 0.5,
                                    true=as.character(ancestry_hat),
                                    false=NA_character_)) %>%
    select(`#FID`, IID, ethnicity=ancestry, ancestry_cluster, ancestry_prob,
           PC1_AVG:PC6_AVG)

ggplot(ancestry_clusters, aes(x=PC1_AVG, PC2_AVG, colour=ancestry_cluster, alpha=ancestry_prob)) +
geom_point()
ggsave('figures/cluster.pc12.png')

ggplot(ancestry_clusters, aes(x=PC3_AVG, PC4_AVG, colour=ancestry_cluster, alpha=ancestry_prob)) +
geom_point()
ggsave('figures/cluster.pc34.png')

ggplot(ancestry_clusters, aes(x=PC5_AVG, PC6_AVG, colour=ancestry_cluster, alpha=ancestry_prob)) +
geom_point()
ggsave('figures/cluster.pc56.png')

write_tsv(ancestry_clusters, 'data/abcd.randomforest.ancestries.tsv')
