#### What is this? ####
# Code written by EN on 2024-02-07 to answer Felix Leiva's request
# The goal is to extract type III marginal sums of squares for fixed and random effects

#### Configuration ####
rm(list=ls())
library(lme4);library(car)

#### Fake data ####
testtemp=rep(c(12,18),25)
gentype=factor(rep(c("A","B"),each=25))
survtime=c(rexp(n = 25,rate = 0.01*testtemp),rexp(n = 25,rate = 0.03*testtemp))
data.fake=data.frame(testtemp,gentype,survtime)
  
#### All data version ####
#Full model
lme0=lmer(survtime~testtemp+(0+gentype+testtemp|gentype),REML = F,data.fake)
summary(lme0)
(ssq.resid <- sum((data.fake$survtime-predict(lme0, re.form=NULL))^2))
ssq.random <- sum((data.fake$survtime - predict(lme0, re.form=NA)) ^ 2)-ssq.resid

# Testemp effect
ssq.fix2 <- sum((predict(lme0, re.form=NA)-mean(predict(lme0, re.form=NA)))^2)

# Marginal gentype effect
lme1=lmer(survtime~testtemp+(0+testtemp|gentype),REML = F,data.fake)
ssq.resid1 <- sum((data.fake$survtime-predict(lme1, re.form=NULL))^2)
(ssq.random1 <- ssq.resid1-ssq.resid)

# Marginal testtemp|gentype effect
lme2=lmer(survtime~testtemp+(0+1|gentype),REML = F,data.fake)
ssq.resid2 <- sum((data.fake$survtime-predict(lme2, re.form=NULL))^2)
(ssq.random2 <- ssq.resid2-ssq.resid)

# Summary full data
ssqs=data.frame(effect=c("testtemp","gentype","testingttem:gentype"),
                typeIII.ssq=c(ssq.fix2,ssq.random1,ssq.random2))
ssqs  

#### Jacknife version #### 
resid.jk=do.call("rbind",lapply(1:nrow(data.fake),function(i)
  {
  #Full model
  lme0=lmer(survtime~testtemp+(0+gentype+testtemp|gentype),REML = F,data.fake[-i,])
  resid.i <- data.fake$survtime[i]-predict(lme0, newdata=data.fake[i,],re.form=NULL)
  random.i <- data.fake$survtime[i] - predict(lme0,newdata=data.fake[i,],re.form=NA)

  # Testemp effect
  fix2.res <- predict(lme0,newdata=data.fake[i,],re.form=NA)-mean(predict(lme0,re.form=NA))

  # Marginal gentype effect
  lme1=lmer(survtime~testtemp+(0+testtemp|gentype),REML = F,data.fake[-i,])
  resid1.i <- data.fake$survtime[i]-predict(lme1, newdata=data.fake[i,],re.form=NULL)
  (random1.i <- resid1.i-resid.i)
  
  # Marginal testtemp|gentype effect
  lme2=lmer(survtime~testtemp+(0+1|gentype),REML = F,data.fake[-i,])
  resid2.i <- data.fake$survtime[i]-predict(lme2, newdata=data.fake[i,], re.form=NULL)
  (random2.i <- resid2.i-resid.i)
  data.frame(fix2.res,random1.i,random2.i)
}))
# Summary jacknife
ssqs.jk=data.frame(effect=c("testtemp","gentype","testingttem:gentype"),
                typeIII.ssq=c(sum(resid.jk$fix2.res^2),sum(resid.jk$random1.i^2),sum(resid.jk$random2.i^2)))
ssqs.jk 


