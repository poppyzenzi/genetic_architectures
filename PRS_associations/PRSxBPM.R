library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(forcats)
library(nnet)
library(gtsummary)
library(magrittr)

## ABCD class data, bpm scores, PRS, associations and plots

##################### load and inspect data ##################### 

# read in bpm scores with class labs - genetic subjects only 
setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/')
bpm_gen_class <- read.csv('./gmm/gmm_abcd/bpm_gen_only_class_labs')

## first check which class is which trajectory
# Convert 'class' and 'eventname' to factors for proper ordering
bpm_gen_class$class <- factor(bpm_gen_class$class)
bpm_gen_class$eventname <- factor(bpm_gen_class$eventname)

# Calculate mean BPM for each 'class' at different 'eventname'
bpm_table <- bpm_gen_class %>%
  group_by(eventname, class) %>%
  summarise(mean_bpm = mean(bpm_sum, na.rm = TRUE)) %>%
  tidyr::pivot_wider(names_from = class, values_from = mean_bpm)

# Print the resulting table
print(bpm_table)

########################################################

#################### make wide df and QC ######################

bpm_gen_class_subset <- subset(bpm_gen_class, select = -c(interview_age))

# make wide
wide_data <- pivot_wider(bpm_gen_class_subset, id_cols = c("src_subject_id", "SubjectNumeric", "sex", "class"),
                         names_from = eventname, values_from = bpm_sum)
# rename cols
colnames(wide_data)[5:ncol(wide_data)] <- paste0("bpm_", colnames(wide_data)[5:ncol(wide_data)])
colnames(wide_data)[1] <- 'IID'

#################### read in and append PRS #################### 
setwd("./prs/prs_bpm_OUT/all_thresholds")
csvs <- c(
  'abcd_adhd_23_prs_0813.t2.best',
  'abcd_asd_prs_0813.t2.best','abcd_bip_prs_0814.t1.best',
  'abcd_mddipsych_prs_0831.t2.5.best',
  #'abcd_mddipsych_prs_sep0907.t2.best',
  'abcd_meta_anx_prs_0813.t1.best',
  'abcd_neu_prs_0813.t1.best', 'abcd_scz_prs_0814.t3.best',
  'bpm_cf_prs_0817.t3.best', 'bpm_high_prs_0817.t2.best',
  'abcd_mood_prs1221.t1.5.best', 'abcd_psychotic_prs1221.t2.best',
  'abcd_neurodev_prs1221.t2.best'
)

# Create an empty list to store dataframes
prs_data_list <- list()

# iterate over each file
for (csv_file in csvs) {
  df <- read.table(csv_file, header = TRUE)
  df <- df[, c('IID', 'PRS')]
  col_prefix <- strsplit(csv_file, "_")[[1]][2]  # extract the prefix of the column name from the file name
  col_prefix2 <- toupper(col_prefix) # capitalise
  setnames(df, c('PRS'), paste0(col_prefix2))  # rename the second column to the appropriate prefix
  prs_data_list[[col_prefix]] <- df  # Store each PRS dataframe in the list
}

# Combine all PRS dataframes into a single dataframe based on IID column
combined_prs_data <- Reduce(function(x, y) merge(x, y, by = 'IID', all.x = TRUE), prs_data_list)
combined_prs_data$ancestry <- 'EUR'

# Merge with smfq_gen_class by IID
merged_data <- merge(wide_data, combined_prs_data, by = 'IID')

# rename PRS columns
setnames(merged_data, c('META', 'CF', 'HIGH','MDDIPSYCH'), c('ANX', 'COMMON', 'HIERARCHICAL','MDD'))

##########################################################

################## principal components #################

PCA <- read.delim('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/multiancestry/abcd.randomforest.ancestries.tsv')
pcs <- PCA[, c('IID', 'PC1_AVG', 'PC2_AVG', 'PC3_AVG', 'PC4_AVG', 'PC5_AVG', 'PC6_AVG')]
df <- merged_data
df <- merge(df, pcs, by = "IID")

####### plot PRS and trajectory association ######

# Using scale() function for Z-score normalization of PRS
for (prs in prs_order) {
  df[[prs]] <- scale(df[[prs]], center = TRUE, scale = TRUE)
}
  
# 1 = persistent, 2 = increasing, 3 = low, 4 = decreasing
# Replace the integers with corresponding labels using recode
df <- df %>% 
  mutate(class = case_when(
    class == 1 ~ "persistent",
    class == 2 ~ "increasing",
    class == 3 ~ "low",
    class == 4 ~ "decreasing",
    TRUE ~ as.character(class) # Add this line if you want to keep other values unchanged
  ))

df$class <- as.factor(df$class)

# select Europeans only
EURdf <- subset(df, ancestry == 'EUR')

# Loop through each PRS variable and fit separate models
individual_models <- list()

for (prs in prs_order) {
  formula <- as.formula(paste("relevel(class, ref = 'low') ~", prs, "+ sex + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG"))
  model <- multinom(formula, data = EURdf)
  
  # Extract coefficients, standard errors, odds ratios, and confidence intervals
  coefficients <- coef(model)
  std_errors <- summary(model)$standard.errors
  odds_ratios <- exp(coefficients)
  lower_ci <- exp(coefficients - 1.96 * std_errors)
  upper_ci <- exp(coefficients + 1.96 * std_errors)
  
  # Store results in a list
  individual_models[[prs]] <- list(
    coefficients = coefficients,
    std_errors = std_errors,
    odds_ratios = odds_ratios,
    upper_ci = upper_ci,
    lower_ci = lower_ci
  )
}

# Create a function to prepare data for plotting
prepare_plot_data <- function(prs_name) {
  odds_ratios <- individual_models[[prs_name]]$odds_ratios
  upper_ci <- individual_models[[prs_name]]$upper_ci
  lower_ci <- individual_models[[prs_name]]$lower_ci
  
  # Get odds ratios for non-reference classes
  plot_data <- data.frame(
    Class = levels(relevel(EURdf$class, ref='low'))[-1],
    Odds_Ratio = odds_ratios,
    PRS = prs_name,
    upper_CI = upper_ci,
    lower_CI = lower_ci
  )

  return(plot_data)
}

## Create an empty list to store the data frames
data_list <- list()

# Loop through each PRS and prepare data for plotting
for (prs in prs_order) {
    prs_data <- prepare_plot_data(prs)
    colnames(prs_data)[colnames(prs_data) == paste0("Odds_Ratio.", prs)] <- "Odds_Ratio.prs"
    colnames(prs_data)[colnames(prs_data) == paste0("upper_CI.", prs)] <- "upper_CI.prs"
    colnames(prs_data)[colnames(prs_data) == paste0("lower_CI.", prs)] <- "lower_CI.prs"
    data_list[[prs]] <- prs_data

}

# Combine the list of data frames into one data frame
combined_data <- do.call(rbind, data_list)

# Plot all PRS on the same plot
# Convert 'Class' to a factor with ordered levels
combined_data$Class <- factor(combined_data$Class, levels = c("decreasing", "increasing", "persistent"))
combined_data$PRS <- factor(combined_data$PRS, levels = c("HIERARCHICAL", "COMMON", "MOOD", "PSYCHOTIC", "NEURODEV", 
                                                          "MDD", "NEU", 
                                                          "ANX", "SCZ", "BIP", 
                                                          "ADHD", "ASD"))

# Define shapes for each level in 'PRS'
shape_vector <- c("HIERARCHICAL" = 17, "COMMON" = 17, "MOOD" = 15, "PSYCHOTIC" = 15, "NEURODEV" = 15, 
                  "MDD" = 16, "NEU" = 16, "ANX" = 16, "SCZ" = 16, "BIP" = 16, "ADHD" = 16, "ASD" = 16)


# Define color vectors for each group
paired_colors <- c("HIERARCHICAL" = "darkgrey", "COMMON" = "grey")
trio1_colors <- c("NEURODEV" = "darkorange3", "PSYCHOTIC" = "deeppink3", "MOOD" = "navy")
trio2_colors <- c("MDD" = "cornflowerblue", "NEU" = "blue", "ANX" = "darkblue")
trio3_colors <- c("SCZ" = "deeppink", "BIP" = "deeppink2")
paired2_colors <- c("ADHD" = "darkorange2","ASD" = "orange")

combined_colors <- c(paired_colors, trio1_colors, trio2_colors, trio3_colors, paired2_colors)

# Plot
ggplot(combined_data, aes(y = Odds_Ratio.prs, x = Class, color = PRS, shape = PRS)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  
  geom_errorbar(aes(ymin = lower_CI.prs , ymax = upper_CI.prs), width = 0.1, position = position_dodge(width = 0.7)) + 
  labs(y = "OR (95% CI)", x = "Trajectory", title = "ABCD EUR") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +  # Set panel background to white
  coord_flip(ylim = c(0.8, 1.7)) +  # Set y-axis limits
  scale_x_discrete() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_shape_manual(values = shape_vector) +
  scale_color_manual(values= combined_colors)

#scale_color_viridis_d(direction = -1)  # Reversed Viridis scale
#scale_color_viridis_d(option = "D")  # Discrete Viridis scale

######################### run the multivariable model ######################### 
## for reviewer 2 (TCE)

multivar_formula <- as.formula(paste("relevel(class, ref = 'low') ~", "MDD + NEU + ANX + SCZ + BIP + ADHD + ASD + sex + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG"))
multivar_model <- multinom(multivar_formula, data = EURdf)

selected_predictors <- c("MDD", "NEU", "ANX", "SCZ", "BIP", "ADHD", "ASD")

## function to make gtsummary::tbl_regression wide format
multinom_pivot_wider <- function(x) {
  # check inputs match expectatations
  if (!inherits(x, "tbl_regression") || !inherits(x$model_obj, "multinom")) {
    stop("`x=` must be class 'tbl_regression' summary of a `nnet::multinom()` model.")
  }
  
  # create tibble of results
  df <- tibble::tibble(outcome_level = unique(x$table_body$groupname_col))
  df$tbl <- 
    purrr::map(
      df$outcome_level,
      function(lvl) {
        gtsummary::modify_table_body(
          x, 
          ~dplyr::filter(.x, .data$groupname_col %in% lvl) %>%
            dplyr::ungroup() %>%
            dplyr::select(-.data$groupname_col)
        )
      }
    )
  
  tbl_merge(df$tbl, tab_spanner = paste0("**", df$outcome_level, "**"))
}

theme_gtsummary_journal(journal = "jama")

## draw the table 
multivar_model %>%
  tbl_regression(exponentiate = T, include=selected_predictors) %>%
  modify_caption('ABCD') %>%
  multinom_pivot_wider  %>%
  modify_table_styling(
    column = c(p.value_1, p.value_2, p.value_3),
    hide = TRUE) %>%
  as_gt() %>%
  gt::tab_options(table.font.names = "Times New Roman")


########################################################################### 


