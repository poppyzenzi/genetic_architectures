library(dplyr)
library(readr)
library(stringr)


demog <- readRDS("/exports/eddie/scratch/s2421111/abcd/acspsw03.rds")
demog <- demog[!duplicated(demog$src_subject_id), ]


ethnic_background <- 
as_tibble(demog) %>%
transmute(FID=rel_family_id, IID=src_subject_id,
          ethnicity=str_replace_all(if_else(is.na(race_ethnicity), true="Do not know", false=as.character(race_ethnicity)),
		                        " ", "_"))


write_tsv(ethnic_background, 'data/ethnic_background.txt')


