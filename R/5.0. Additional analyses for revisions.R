# ------------------------------------------------------------------------------
# This script estimate CTmax and z for after fitting a linear model of survival
# time and test temperature for each DGRP line and sex as well as their errors
# and CI for both estimates, for each genotype and sex. I hope this code will be
# useful to other as well working in similar subjects
# ------------------------------------------------------------------------------
# Script created by Felix P. Leiva
# Created on: 20241107
# Modifications: XXXXX by YYYYYYY

# If you use it, please cite as:
# Leiva FP, Santos M, Niklitschek EJ, Rezende EL, & Verberk WCEP. (2024). Paper
# data and code for: Genetic variation of heat tolerance in a model ectotherm:
# an approach using thermal death time curves. Zenodo.
# https://doi.org/10.5281/zenodo.12155988.

# Cleaning working space
rm(list=ls())

set.seed(6955)

# check directory
getwd()
# ------------------------------------------------------------------------------
#Libraries
library(dplyr)
library(MASS)
library(broom)
library(ggplot2)
library(msm)  # For deltamethod function
library(boot) # For bootstrap
# ------------------------------------------------------------------------------
# Load data of the 20 DGRP lines ####
all.lines.21<-read.csv("../Outputs/0.1.1. Survival time of Drosophila melanogaster across DGRP lines.csv")
# check data structure
str(all.lines.21)
#stock as factor
all.lines.21$stock<-as.factor(all.lines.21$stock)
#transform minutes in log10-scale to minutes
all.lines.21$surv.time2<-as.numeric(10^(all.lines.21$surv.time))
# subset by sex
male <- subset (all.lines.21, all.lines.21$sex == "male")
female <- subset (all.lines.21, all.lines.21$sex != "male")
# ------------------------------------------------------------------------------
# I just copied below the reviewer comment to understand better what he/she wants

# I thought it would be helpful to report the variance associated with the
# estimates of CTmax and z, especially because the authors are using linear
# regression that extrapolates beyond the range of measurements. The authors
# calculated the standard error across all lines measured, but I think it would
# be more informative to report the standard error of each line (and if possible,
# incorporate this uncertainty in downstream analyses). The model output should
# report the standard errors associated with B0 and B1, and there should also be
# a way to calculate the standard error associated with an estimated value
# derived from the model (i.e., the CTmax estimate). Reporting these standard
# errors, and/or confidence intervals for each estimate, would help indicate the
# precision of these estimates within each line.

# To address the reviewer's comment, we can calculate and report the standard
# errors and confidence intervals for CTmax and z for each DGRP line. Here's how
# we can modify the analysis to include this information:
  
# First, let's create a function to fit the model and extract coefficients along with their standard errors:
TDT_curve_with_se <- function(df) {
  model <- lm(surv.time ~ test.temp, data = df)
  coef_summary <- summary(model)$coefficients
  
  b0 <- coef_summary[1, 1]    # extract intercept
  b1 <- coef_summary[2, 1]    # extract slope
  se_b0 <- coef_summary[1, 2] # extract SE intercept
  se_b1 <- coef_summary[2, 2] # extract SE slope
  
  # # Calculate CTmax and its standard error using the delta method. More of the delta methods can be found here:
  # https://bookdown.org/ts_robinson1994/10EconometricTheorems/dm.html#applied-example
  ctmax <- -b0 / b1
  se_ctmax <- deltamethod(~ -x1/x2, c(b0, b1), vcov(model))
  
  # Calculate z and its standard error
  z <- -1 / b1
  se_z <- deltamethod(~ -1/x2, c(b0, b1), vcov(model))
  
  data.frame(
    b0 = round(b0, 4),
    se_b0 = round(se_b0, 4),
    b1 = round(b1, 4),
    se_b1 = round(se_b1, 4),
    CTmax = round(ctmax, 4),
    se_CTmax = round(se_ctmax, 4),
    z = round(z, 4),
    se_z = round(se_z, 4)
  )
}

#  Now lets apply this function to our data:
  
# for females
results_female_OLR <- female %>%
  group_by(genotype) %>%
  do(TDT_curve_with_se(.)) %>%
  ungroup()

# for males
results_male_OLR <- male %>%
  group_by(genotype) %>%
  do(TDT_curve_with_se(.)) %>%
  ungroup()

# Calculate 95% confidence intervals:

calculate_ci <- function(estimate, se) {
  lower <- estimate - 1.96 * se
  upper <- estimate + 1.96 * se
  paste0("[", round(lower, 4), ", ", round(upper, 4), "]")
}

results_female_OLR <- results_female_OLR %>%
  mutate(
    CI_CTmax = calculate_ci(CTmax, se_CTmax),
    CI_z = calculate_ci(z, se_z)
  )

results_male_OLR <- results_male_OLR %>%
  mutate(
    CI_CTmax = calculate_ci(CTmax, se_CTmax),
    CI_z = calculate_ci(z, se_z)
  )

# Create and export Table S2:

names(results_female_OLR) <- c("DGRP line", 
                               "b0_female", "se_b0_female", "b1_female", "se_b1_female", 
                               "CTmax_female", "se_CTmax_female", "z_female", "se_z_female",
                               "CI_CTmax_female", "CI_z_female")
names(results_male_OLR) <- c("DGRP line_male", 
                             "b0_male", "se_b0_male", "b1_male", "se_b1_male", 
                             "CTmax_male", "se_CTmax_male", "z_male", "se_z_male",
                             "CI_CTmax_male", "CI_z_male")

TableS2 <- cbind(results_female_OLR, results_male_OLR) 
names(TableS2)

select_columns <- c("DGRP line",
                    "b0_female", "se_b0_female", "b1_female", "se_b1_female",
                    "CTmax_female", "se_CTmax_female", "z_female", "se_z_female",
                    "CI_CTmax_female", "CI_z_female",
                    "b0_male", "se_b0_male", "b1_male", "se_b1_male",
                    "CTmax_male", "se_CTmax_male", "z_male", "se_z_male",
                    "CI_CTmax_male", "CI_z_male"
)

# Seleccionar las columnas
TableS2 <- TableS2[, select_columns]

write.csv(TableS2, "../Outputs/Table_S2_CTmax_and_z_per_DGRP_line_and_sex_with_SE_and_CI.csv", row.names = FALSE)

# ... (previous code remains the same)

# Bootstrap function
bootstrap_TDT <- function(data, indices) {
  d <- data[indices,]
  model <- lm(surv.time ~ test.temp, data = d)
  coef <- coef(model)
  b0 <- coef[1]
  b1 <- coef[2]
  ctmax <- -b0 / b1
  z <- -1 / b1
  c(ctmax = ctmax, z = z)
}

# Function to perform bootstrap and calculate CI and SE
bootstrap_analysis <- function(df) {
  set.seed(123)  # for reproducibility
  boot_results <- boot(data = df, statistic = bootstrap_TDT, R = 1000)
  
  ctmax_ci <- boot.ci(boot_results, type = "perc", index = 1)
  z_ci <- boot.ci(boot_results, type = "perc", index = 2)
  
  data.frame(
    CTmax = mean(boot_results$t[,1]),
    CTmax_SE = sd(boot_results$t[,1]),
    CTmax_lower_CI = ctmax_ci$percent[4],
    CTmax_upper_CI = ctmax_ci$percent[5],
    z = mean(boot_results$t[,2]),
    z_SE = sd(boot_results$t[,2]),
    z_lower_CI = z_ci$percent[4],
    z_upper_CI = z_ci$percent[5]
  )
}

# Apply bootstrap analysis to each genotype for females and males
results_female_bootstrap <- female %>%
  group_by(genotype) %>%
  do(bootstrap_analysis(.)) %>%
  ungroup()

results_male_bootstrap <- male %>%
  group_by(genotype) %>%
  do(bootstrap_analysis(.)) %>%
  ungroup()

# Combine results
names(results_female_bootstrap) <- c("DGRP line", 
                                     "CTmax_female", "CTmax_SE_female", "CTmax_lower_CI_female", "CTmax_upper_CI_female",
                                     "z_female", "z_SE_female", "z_lower_CI_female", "z_upper_CI_female")
names(results_male_bootstrap) <- c("DGRP line", 
                                   "CTmax_male", "CTmax_SE_male", "CTmax_lower_CI_male", "CTmax_upper_CI_male",
                                   "z_male", "z_SE_male", "z_lower_CI_male", "z_upper_CI_male")

TableS2_bootstrap <- full_join(results_female_bootstrap, results_male_bootstrap, by = "DGRP line")

# Round numeric columns to 4 decimal places
TableS2_bootstrap <- TableS2_bootstrap %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Export results
write.csv(TableS2_bootstrap, "../Outputs/Table_S2_CTmax_and_z_per_DGRP_line_and_sex_bootstrap_with_SE.csv", row.names = FALSE)

# Optional: Create a summary table with mean, SE, and CI across all lines
summary_stats <- TableS2_bootstrap %>%
  summarise(across(where(is.numeric), list(mean = mean, se = ~sd(.)/sqrt(n()))))

write.csv(summary_stats, "../Outputs/Summary_stats_CTmax_and_z_across_lines.csv", row.names = FALSE)