# ------------------------------------------------------------------------------
# This script adapt a matlab to R. Code estimate heritability for heat tolerance
# of DGRP lines. 

# Originally, the matlab code was written by Mauro Santos (UAB,
# Spain)

# This R script was created by Edwin j. Niklitschek

# Cleaning working space
rm(list=ls())


set.seed(42)

# check directory
getwd()
# ------------------------------------------------------------------------------
#Libraries
library(car)
library(lme4)
library(MuMIn)
library(dplyr)
library(nlme)
# +++++++++++++++++++++++++++++++++++++++
#### Load data of the 20 DGRP lines ####
all.lines.21<-read.csv("../Outputs/0.1.1. Survival time of Drosophila melanogaster across DGRP lines.csv")
# check data structure
str(all.lines.21)
#stock as factor
all.lines.21$stock<-as.factor(all.lines.21$stock)
#transform minutes in log10-scale to minutes
all.lines.21$surv.time2<-as.numeric(10^(all.lines.21$surv.time))
# subset by sex
names(all.lines.21)[c(14,5,2)]=c("Survtime","Testtemp","Genotype")
Data_females<-subset(all.lines.21,sex=="female") #females
Data_males<-subset(all.lines.21,sex=="male") #males
n_0Females <- 10.2479840885613
n_0Males = 10.8054730372430

#### Ancova ####
ancova=function(data,n_0){
  ancova1=lm(Survtime ~ -1+Genotype/Testtemp, data = data)
  # Extract coefficients
  coefficients <- summary(ancova1)$coefficients
  # Intercepts and slopes
  intercepts <- coefficients[1:length(unique(data$Genotype)), "Estimate"]
  slopes <- coefficients[(length(unique(data$Genotype))+1):(length(unique(data$Genotype))*2), "Estimate"]
  Mean_b0 = mean(intercepts); Mean_b1 = mean(slopes);
  # Calculos
  vr=matrix(c(var(intercepts),rep(cov(intercepts,slopes),2),var(slopes)),nr=2)
  # Replicando a Mauro
  Mean_b0F = mean(intercepts); Mean_b1F = mean(slopes);
  Var_b0F = as.numeric(vr[1,1]); Var_b1F = as.numeric(vr[2,2]);
  Cov_b0F_b1F = as.numeric(vr[2,1])
  Var_eF = sigma(ancova1)^2;
  
  # Calculate variance components
  Var_CTmaxF <- ((Mean_b0F)^2 / (Mean_b1F)^2) * 
    (Var_b0F/(Mean_b0F)^2 - (2 * Cov_b0F_b1F / (Mean_b0F * Mean_b1F)) + Var_b1F/(Mean_b1F)^2)
  Var_zF <- Var_b1F / (Mean_b1F)^4
  
  # Calculate heritability estimates
  n_0
  H2_CTmaxF <- Var_CTmaxF / (Var_CTmaxF + Var_zF + Var_eF / n_0)
  H2_zF <- Var_zF / (Var_CTmaxF + Var_zF + Var_eF / n_0)
  as.data.frame(t(data.frame(Var_b0F, Var_b1F, Cov_b0F_b1F, Var_eF, Var_CTmaxF, Var_zF, H2_CTmaxF, H2_zF)))
}  

#### Funcion LME ####
lme2=function(data,n_0,method){
  #1. All data estimates
  lme.f <- lme(Survtime ~ Testtemp,random=list(~Testtemp|Genotype,~1|Genotype), data = data, method=method,control = lmeControl(maxiter=1000) )
  (vr=VarCorr(lme.f))
  # Replicando a Mauro
  Mean_b0F = fixef(lme.f)[1]; Mean_b1F = fixef(lme.f)[2];
  Var_b0F = as.numeric(vr[2,1]); Var_b1F = as.numeric(vr[3,1]);
  Cov_b0F_b1F = as.numeric(vr[2,2])*as.numeric(vr[3,2])*as.numeric(vr[3,3])
  Var_eF = sigma(lme.f)^2;
  
  # Calculate variance components
  Var_CTmaxF <- ((Mean_b0F)^2 / (Mean_b1F)^2) * 
    (Var_b0F/(Mean_b0F)^2 - (2 * Cov_b0F_b1F / (Mean_b0F * Mean_b1F)) + Var_b1F/(Mean_b1F)^2)
  Var_zF <- Var_b1F / (Mean_b1F)^4
  
  # Calculate heritability estimates
  n_0
  H2_CTmaxF <- Var_CTmaxF / (Var_CTmaxF + Var_zF + Var_eF / n_0)
  H2_zF <- Var_zF / (Var_CTmaxF + Var_zF + Var_eF / n_0)
  as.data.frame(t(data.frame(Var_b0F, Var_b1F, Cov_b0F_b1F, Var_eF, Var_CTmaxF, Var_zF, H2_CTmaxF, H2_zF)))
}
### Funcion Jacknife ####
jk=function(data,n_0,method,funcion){
    if(funcion=="ancova"){
      jk.all=ancova(data = data,n_0=n_0)
    } else {if(funcion=="lme") jk.all=lme2(data=data,n_0=n_0,method=method) else {
      print('Incorrect function specification. Must be either "ancova" or "lme"')}}
  jk.i=do.call("cbind",lapply(unique(data$Genotype), function(i){
    data.i=subset(data,Genotype!=i)
    if(funcion=="ancova"){
      ancova(data = data.i,n_0=n_0)
    } else {if(funcion=="lme") lme2(data=data.i,n_0=n_0,method=method) else {
      print('Incorrect function specification. Must be either "ancova" or "lme"')}}
  }))
  ndata=length(unique(data$Genotype))
  pseudo.values=do.call("cbind",apply(jk.i,2,function(i) ndata*jk.all-(ndata-1)*i))
  Jack_estimates=matrix(rowMeans(pseudo.values),nc=1)
  #dv.i=do.call("cbind",apply(pseudo.values,1, function(i) i-jk.all))
  dv.i=apply(pseudo.values,2, function(i) i-Jack_estimates)
  SE=sqrt(rowSums(dv.i^2)/(ndata*(ndata-1)))
  low = Jack_estimates - 2*SE;
  up = Jack_estimates + 2*SE;
  jk=(data.frame(Jacknife=Jack_estimates,low,up))
}

#### Aplica funciones hembras ####
#Hembras Ancova
females.ancova=ancova(data = Data_females,n_0=n_0Females)
females.jk.ancova=jk(data=Data_females,n_0=n_0Females,funcion = "ancova")
females.ancova=cbind(females.ancova,females.jk.ancova)
#Hembras ML
females.est.ML=lme2(data = Data_females,n_0=n_0Females,method="ML")
females.jk.ML=jk(data=Data_females,n_0=n_0Females,method="ML",funcion = "lme")
females.ML=cbind(females.est.ML,females.jk.ML)
#Hembras REML
females.est.REML=lme2(data = Data_females,n_0=n_0Females,method="REML")
females.jk.REML=jk(data=Data_females,n_0=n_0Females,method="REML",funcion="lme")
females.REML=cbind(females.est.REML,females.jk.REML)

#### Aplica funciones machos ####
#Machos Ancova
males.ancova=ancova(data = Data_males,n_0=n_0Males)
males.jk.ancova=jk(data=Data_males,n_0=n_0Males,funcion = "ancova")
males.ancova=cbind(males.ancova,males.jk.ancova)
#Machos ML
males.est.ML=lme2(data = Data_males,n_0=n_0Males,method="ML")
males.jk.ML=jk(data=Data_males,n_0=n_0Males,method="ML",funcion="lme")
males.ML=cbind(males.est.ML,males.jk.ML)
#Machos REML
males.est.REML=lme2(data = Data_males,n_0=n_0Males,method="REML")
males.jk.REML=jk(data=Data_males,n_0=n_0Males,method="REML",funcion="lme")
males.REML=cbind(males.est.REML,males.jk.REML)

# combine
ancova <- cbind(females.ancova, males.ancova)
ML     <- cbind(females.ML, males.ML)
REML   <- cbind(females.REML, males.REML)

names(ancova)
colnames(ancova) <- c("Female_estimate", "Females_Jacknife", "Females_low_95%_CI", "Females_up_95%_CI", 
                      "Males_estimate", "Males_Jacknife", "Males_low_95%_CI", "Males_up_95%_CI")
ancova$Method <- "Anova"

names(ML)
colnames(ML) <- c("Female_estimate", "Females_Jacknife", "Females_low_95%_CI", "Females_up_95%_CI", 
                      "Males_estimate", "Males_Jacknife", "Males_low_95%_CI", "Males_up_95%_CI")
ML$Method <- "ML"

names(REML)
colnames(REML) <- c("Female_estimate", "Females_Jacknife", "Females_low_95%_CI", "Females_up_95%_CI", 
                  "Males_estimate", "Males_Jacknife", "Males_low_95%_CI", "Males_up_95%_CI")
REML$Method <- "REML"

# Combine and export to csv
Table_1 <- rbind(ancova, ML, REML)
write.csv(Table_1, "Estimates of variance-covariance components and broad-sense heritability DGRP lines.csv", row.names = TRUE)
