library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(forcats)
library(nnet)
library(gtsummary)
library(magrittr)

## ABCD class data, bpm scores, PRS: sensitivity and attrition analyses

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

############################################################### 

## regression of each PRS with smfq_1-4

df <- merged_data

# Initialize an empty dataframe to store the regression results
results <- data.frame()

# Loop through each PRS column
for (i in 13:21) {
  prs_name <- colnames(df)[i]
  # Loop through each bpm column
  for (j in 5:12) {
    bpm_name <- colnames(df)[j]
    
    # Run linear regression
    lm_result <- lm(scale(df[, j]) ~ scale(df[, i]))  # Standardize predictors and outcome
    
    # Store the results in a dataframe
    result_row <- data.frame(
      PRS = prs_name,
      BPM_Time = bpm_name,
      Coefficients = coef(lm_result)[2],
      p_value = summary(lm_result)$coefficients[2, 4]
    )
    # Append the results to the overall dataframe
    results <- rbind(results, result_row)
  }
}


# Calculate confidence intervals for coefficients
results$lower_ci <- NA
results$upper_ci <- NA

for (index in 1:nrow(results)) {
  prs <- results$PRS[index]
  bpm <- results$BPM_Time[index]
  # Subset data for each PRS and SMFQ combination
  subset_data <- df[, c(prs, bpm)]
  subset_data <- subset_data[complete.cases(subset_data), ]  # Remove rows with NA
  # Run linear regression
  lm_result <- lm(scale(subset_data[, 2]) ~ scale(subset_data[, 1]))  # Standardize predictors and outcome
  # Calculate confidence intervals
  conf_int <- confint(lm_result)[2, ]
  results$lower_ci[index] <- conf_int[1]
  results$upper_ci[index] <- conf_int[2]
}

# Print the table of regression results
print(results)


######## Plot the results ###########

# Specify the order of the PRS variable
prs_order <- c("ANX", "NEU", "MDD", "BIP", "SCZ", "ADHD", "ASD", "MOOD", "PSYCHOTIC", "NEURODEV", "COMMON", "HIERARCHICAL")

# Convert PRS to a factor with the desired order
results$PRS <- factor(results$PRS, levels = prs_order)

# plot with confidence intervals
ggplot(results, aes(x = PRS, y = Coefficients, color = BPM_Time)) +
  geom_point(position = position_dodge(width = 0.5)) +  # Adding separation between points
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.5)) +  # Separation between error bars
  labs(title = "Coefficients of PRS vs. BPM_Time",
       x = "Polygenic Risk Score (PRS)",
       y = "Betas") +
  theme_minimal()

######## now do the same for number of responses #########

# Calculate the sum of non-NA values in columns 5 to 8 and create a new column 'responses'
merged_data$missed_surveys <- rowSums(is.na(merged_data[, 5:12]))

# Initialize an empty dataframe to store the regression results
df2 <- merged_data
# Initialize an empty dataframe to store the regression results
results2 <- data.frame()

# Loop through each PRS column
for (i in 13:24) {
  prs_name <- colnames(df2)[i]
  # Subset the data to rows with non-NA responses for the current PRS and response columns
  subset_data <- df2[!is.na(df2$missed_surveys) & !is.na(df2[, i]), ]
  # Standardize predictors and response variable
  scaled_data <- scale(subset_data[, c('missed_surveys', prs_name)])
  # Run linear regression
  lm_result <- lm(scaled_data[, 1] ~ scaled_data[, 2])
  # Store the results in a dataframe
  result_row <- data.frame(
    PRS = prs_name,
    Beta_Coefficient = coef(lm_result)[2],
    p_value = summary(lm_result)$coefficients[2, 4],
    # Calculate confidence intervals
    lower_ci = confint(lm_result)[2, 1],
    upper_ci = confint(lm_result)[2, 2]
  )
  # Append the results to the overall dataframe
  results2 <- rbind(results2, result_row)
}

# Print the table of regression results
print(results2)

results2$PRS <- factor(results2$PRS, levels = prs_order)

# plot results2 with error bars
ggplot(results2, aes(x = PRS, y = Beta_Coefficient)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  labs(title = "PRS vs. # missed surveys", x = "Polygenic Risk Score (PRS)", y = "Beta Coefficient")

########## inspect missing data (attrition) ###############

# Count rows with 2 or more missing values in smfq_1 to smfq_4 columns

missing_time_points <- data.frame()

for (i in 1:4) {
  missing_n_or_more <- rowSums(is.na(merged_data[, 5:12])) >= i
  x <- sum(missing_n_or_more)
  missing_time_points <- rbind(missing_time_points, data.frame(TimePoints = i, Count = x))
}

missing_time_points


##########################################################################
################### PRS x symptom level data ############################


setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
bpm_symptoms <- read.table('abcd_bpm_sym_long.txt')
colnames(bpm_symptoms)[colnames(bpm_symptoms) == "src_subject_id"] <- "IID"

# merge combined PRS data 
symp_with_PRS <- merge(bpm_symptoms, combined_prs_data, by='IID')

setnames(symp_with_PRS, c('MDDIPSYCH', 'META', 'CF', 'HIGH'), c('MDD', 'ANX', 'COMMON', 'HIERARCHICAL'))

setnames(symp_with_PRS, c('bpm_9_y','bpm_11_y','bpm_12_y','bpm_13_y','bpm_18_y','bpm_19_y'),
         c('worthless or inferior','fearful or anxious', 'feel guilty', 'self-conscious',
           'unhappy, sad, depressed', 'worry'))

# symptoms are already binary

symptom_columns <- 3:8
prs_columns <- 9:15

# Initialize an empty dataframe to store the regression results
results <- data.frame()

for (i in 9:15) {
  prs_name <- colnames(symp_with_PRS)[i]
  for (j in 3:8) {
    bpm_name <- colnames(symp_with_PRS)[j]
    lm_result <- lm((symp_with_PRS[, j]) ~ scale(symp_with_PRS[, i]))  # Standardize predictors and outcome
    
    # Store the results in a dataframe
    result_row <- data.frame(
      PRS = prs_name,
      BPM_Symptom = bpm_name,
      Coefficients = coef(lm_result)[2],
      p_value = summary(lm_result)$coefficients[2, 4]
    )
    # Append the results to the overall dataframe
    results <- rbind(results, result_row)
  }
}

# Calculate confidence intervals for coefficients
results$lower_ci <- NA
results$upper_ci <- NA

for (index in 1:nrow(results)) {
  prs <- results$PRS[index]
  bpm <- results$BPM_Symptom[index]
  # Subset data for each PRS and SMFQ combination
  subset_data <- symp_with_PRS[, c(prs, bpm)]
  subset_data <- subset_data[complete.cases(subset_data), ]  # Remove rows with NA
  # Run linear regression
  lm_result <- lm((subset_data[, 2]) ~ scale(subset_data[, 1]))  # Standardize predictors and outcome
  # Calculate confidence intervals
  conf_int <- confint(lm_result)[2, ]
  results$lower_ci[index] <- conf_int[1]
  results$upper_ci[index] <- conf_int[2]
}

# Specify the order of the PRS variable
prs_order <- c("ANX", "NEU", "MDD", "BIP", "SCZ", "ADHD", "ASD"
               #, "MOOD", "PSYCHOTIC", "NEURODEV", "COMMON", "HIERARCHICAL"
               )
symptom_order <-  c('worthless or inferior','fearful or anxious', 'feel guilty', 'self-conscious',
                    'unhappy, sad, depressed', 'worry')


# Convert to a factor with the desired order
results$PRS <- factor(results$PRS, levels = prs_order)
results$BPM_Symptom <- factor(results$BPM_Symptom, levels = symptom_order)

theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggplot(results, aes(x = Coefficients, y = BPM_Symptom, color = PRS)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.5)) +
  labs(title = "Coefficients of PRS vs. BPM_Symptom",
       x = "Betas",
       y = "BPM Symptom") +
  theme_minimal() 
