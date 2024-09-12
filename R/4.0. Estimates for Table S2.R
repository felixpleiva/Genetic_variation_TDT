# ------------------------------------------------------------------------------
# This script estimate CTMAX and z for each DGRP line and sex

# Cite as:

# Leiva FP, Santos M, Niklitschek EJ, Rezende EL, & Verberk WCEP. (2024). Paper
# data and code for: Genetic variation of heat tolerance in a model ectotherm:
# an approach using thermal death time curves. Zenodo.
# https://doi.org/10.5281/zenodo.12155988.

# Cleaning working space
rm(list=ls())

set.seed(42)

# check directory
getwd()
# ------------------------------------------------------------------------------
#Libraries
library(dplyr)
library(MASS)
# ------------------------------------------------------------------------------
#### Load data of the 20 DGRP lines ####
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

#-------------------------------------------------------------------------------
# Function to fit standard linear regression and extract coefficients (b0 = intercept; b1: slope)
fit_model <- function(df) {
  model <- lm(surv.time ~ test.temp, data = df)
  coef <- coef(model)
  data.frame(b0 = round(coef[1], 4), b1 = round(coef[2], 4))
  }

# Group data by genotype, fit model, and extract coefficients

results_female_OLR <- female %>%
  group_by(genotype) %>%
  do(fit_model(.)) %>%
  ungroup()

results_male_OLR <- male %>%
  group_by(genotype) %>%
  do(fit_model(.)) %>%
  ungroup()

# calculate CTmax
results_female_OLR$CTmax <- round((results_female_OLR$b0/results_female_OLR$b1) * -1, 4)
mean(results_female_OLR$CTmax) # 41.99027
sd(results_female_OLR$CTmax)   # 0.7709354

results_male_OLR$CTmax   <- round((results_male_OLR$b0/results_male_OLR$b1) * -1, 4)
mean(results_male_OLR$CTmax) # 43.19402
sd(results_male_OLR$CTmax)   # 1.809911

# calculate Z
results_female_OLR$z <- round(-1/results_female_OLR$b1, 4)
mean(results_female_OLR$z)   # 2.75906
sd(results_female_OLR$z)     # 0.4269119

results_male_OLR$z   <- round(-1/results_male_OLR$b1, 4)
mean(results_male_OLR$z)   # 3.559465
sd(results_male_OLR$z)     # 1.122591
#-------------------------------------------------------------------------------
# Function to fit robust linear regression and extract coefficients (b0 = intercept; b1: slope)
fit_model <- function(df) {
  model <- rlm(surv.time ~ test.temp, 
               psi = MASS::psi.bisquare,
               c = 4.685,
               scale.est = "MAD",
               data = df)
  coef <- coef(model)
    data.frame(
    b0 = round(coef[1], 4),
    b1 = round(coef[2], 4)
  )
}

# Group data by genotype, fit model, and extract coefficients

results_female_RLR <- female %>%
  group_by(genotype) %>%
  do(fit_model(.)) %>%
  ungroup()

results_male_RLR <- male %>%
  group_by(genotype) %>%
  do(fit_model(.)) %>%
  ungroup()

# calculate CTmax
results_female_RLR$CTmax <- round((results_female_RLR$b0/results_female_RLR$b1) * -1, 4)
results_male_RLR$CTmax   <- round((results_male_RLR$b0/results_male_RLR$b1) * -1, 4)

# calculate Z
results_female_RLR$z <- round(-1/results_female_RLR$b1, 4)
results_male_RLR$z   <- round(-1/results_male_RLR$b1, 4)

#-------------------------------------------------------------------------------
# Export Table S2 (only result from OLS)
TableS2 <- cbind(results_female_OLR, results_male_OLR)
names(TableS2)
TableS2 = subset(TableS2, select = -(6))
names(TableS2)
names(TableS2) <- c("DGRP line", 
                    "b0_female", "b1_female", "CTmax_female", "z_female",
                    "b0_male", "b1_male", "CTmax_male", "z_male")

TableS2
# export Table S2
write.csv(TableS2, "../Outputs/Table_2 CTmax and z per DGRP line and sex.csv", row.names = F)    
