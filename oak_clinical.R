library(rjags)
library(MASS)
library(broom)
library(tidyverse)
library(ggplot2)
library(parallel)
library(R2jags)
library(bayesplot)
library(dplyr)
library(runjags)

load('.RData')

options(scipen = 999)

####################################################################
## introducing the dataset
####################################################################
data = read.csv('oak_clinical.csv', header = T)
str(data)
head(data)


##in total we will have five covariates: age, sex, smoking, body condition, TC/IC score
##which correspond to BAGE, SEX, TOBHX, ECOGGR, and TC1IC1.

##the treatment variable is TRT01P. We set "Docetaxel" = 0 and "MPDL3280A" = 1

##the response variable is OS, and censoring indicator OS.CNSR

##Set up the X matrix
X = cbind(data$TRT01P, data$BAGE, data$SEX, data$TOBHX, data$ECOGGR, data$TC1IC1)
X = as.data.frame(X)
colnames(X) = c('trt', 'age', 'sex', 'smoke', 'body', 'tcic')
##response
y = cbind(data$OS, data$OS.CNSR)
colnames(y) = c('survT', 'censor')

############
##Cleaning the dataset
############

##treatment variable
X[,1] = ifelse(X[,1] == 'Docetaxel', 0, 1)

##Transform it to integers
X[,2] = as.numeric(X[,2])
X[,3] = ifelse(X[,3] == 'M', 0, 1)
X[,5] = as.numeric(X[,5])

##age, sex, smoking and body condition are all characters. No missing occurs.

##TC/IC score
tcic.index = unique(X$tcic)
for(i in 1:length(X$tcic)){
  if(X[i, 6] == tcic.index[1]){
    X[i, 6] = 0
  }else if(X[i, 6] == tcic.index[2]) {
    X[i, 6] = 1
  }else{
    X[i, 6] = 2}
}
X[,6] = as.numeric(X[,6])

##change smoke into non-smoke = 0, smoked = 1
for(i in 1:length(X$smoke)){
  if(X[i, 4] == 'NEVER'){
    X[i,4] = 0
  }else{
    X[i,4] = 1
  }
}
X[,4] = as.numeric(X[,4])


##remove UNKNOWNS
df = cbind(y, X)
df = df[!grepl(2, df$tcic), ] ##remove UNKNOWN

##normalize age
df$age = (df$age - mean(df$age))/sd(df$age)





####################################################################
## Visualization
####################################################################
library(survival)
library(ggplot2)
library(GGally)

df$dead = df$censor
##change from censoring indicator to a death indicator
for(i in 1:nrow(df)){
  if(df$censor[i] == 1){
    df$dead[i] = 0
  }else{
    df$dead[i] = 1
  }
}

df$Trt = I(df$trt)
sf = survfit(Surv(survT, dead) ~ Trt, data=df)
ggsurv(sf) + geom_vline(xintercept = 2.7) + coord_cartesian(xlim = c(0, 30)) + 
  theme(plot.title = element_text(hjust = 0.6, face = 'bold', size = 20))


####################################################################
## Fitting the models
####################################################################

##First, we are going to fit a standard Bayesian AFT model, assuming Weibull

modelstring.weibull = '
model{
  for(i in 1 : N) {
    is.censored[i] ~ dinterval(time[i], centime[i])
    time[i] ~ dweib(b, lambda[i])
    lambda[i] <- exp(-(mu + beta_formula[i]) * b)
    beta_formula[i] <- beta[1] * trt[i] + 
                       beta[2] * age[i] + 
                       beta[3] * sex[i] +
                       beta[4] * body[i] +
                       beta[5] * tcic[i] + 
                       beta[6] * smoke[i]
  }
  
  ##priors for betas
  mu ~ dnorm(0, 0.001)
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for b
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.weibull = file.path(dir, 'weibull.txt')
writeLines(modelstring.weibull, con=modelfile.weibull)

##we need to modify survT for using dweibull in JAGS. Separate it into mod.time and
##mod.centime
df$mod.time = df$survT
df$centime = df$survT
for(i in 1:nrow(df)){
  if(df$censor[i] == 0) {
    ##if no censoring 
    df$centime[i] = 10000
  }else{
    ##if censoring
    df$mod.time[i] = NA
  }
}

datlist = list(N = nrow(df),
               time = df$mod.time,
               centime = df$centime,
               is.censored = df$censor,
               trt = df$trt,
               age = df$age,
               sex = df$sex,
               body = df$body,
               tcic = df$tcic,
               smoke = df$smoke)

weibull.prior = jags.parallel(data = datlist,
                              parameters.to.save = c('beta', 'mu', 'sigma'),
                              n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                              n.thin = 1, model.file = modelfile.weibull)

weibull.prior$BUGSoutput[10]$summary[,c(1,3,7)]
#plot(weibull.prior)
#traceplot(weibull.prior)


modelstring.weibull.cus = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
    
    beta_formula[i] <- beta[1] * trt[i] + 
                       beta[2] * age[i] + 
                       beta[3] * sex[i] +
                       beta[4] * body[i] +
                       beta[5] * tcic[i] + 
                       beta[6] * smoke[i]
    
    h[i] <- (b/time[i]) * exp(b * (log(time[i]) - mu - beta_formula[i]))
    #lambda[i] * b * pow(time[i], b-1)
    S[i] <- exp(-exp(b * (log(time[i]) - mu - beta_formula[i])))
    #exp(-(lambda[i]) * pow(time[i], b))
    lambda[i] <- exp(-(mu + beta_formula[i]) * b)
  }
  
  ##priors for betas
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.weibull.cus = file.path(dir, 'weibullcus.txt')
writeLines(modelstring.weibull.cus, con=modelfile.weibull.cus)

datlist.tau = list(N = nrow(df),
                   time = df$survT,
                   is.observed = df$dead,
                   trt = df$trt,
                   age = df$age,
                   sex = df$sex,
                   body = df$body,
                   tcic = df$tcic,
                   smoke = df$smoke,
                   pi = pi)

weibull.prior.cus = jags.parallel(data = datlist.tau, 
                                  parameters.to.save = c('beta', 'mu', 'sigma'),
                                  n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                                  n.thin = 1, model.file = modelfile.weibull.cus)

weibull.prior$BUGSoutput[10]$summary[,c(1,3,7)]
weibull.prior.cus$BUGSoutput[10]$summary[,c(1,3,7)]


library(SurvRegCensCov)
freq.weibull = WeibullReg(formula = Surv(survT, dead) ~ I(trt) + age + I(sex) +
                            I(body) + I(tcic) + I(smoke.prev) + I(smoke.now), data = df)
freq.weibull
weibull.prior$BUGSoutput[10]$summary[,c(1,3,7)]
weibull.prior.cus$BUGSoutput[10]$summary[,c(1,3,7)]




modelstring.logn = '
model{
  for(i in 1 : N) {
    is.censored[i] ~ dinterval(time[i], centime[i])
    time[i] ~ dlnorm(mean[i], b)
    mean[i] <- mu + beta[1] * trt[i]
  }
  
  ##priors for betas
  beta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1/2)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.logn = file.path(dir, 'logn.txt')
writeLines(modelstring.logn, con=modelfile.logn)

ln.prior = jags.parallel(data = datlist,
                         parameters.to.save = c('beta', 'mu', 'sigma'),
                         n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                         n.thin = 1, model.file = modelfile.logn)
ln.prior$BUGSoutput[10]$summary[,c(1,3,7)]
traceplot(ln.prior)


modelstring.logn.cus = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
                       h.norm[i]*S.norm[i],
                       S.norm[i]))
    
    ##constructing hazard and Survival function
    
    epsilon.form.norm[i] <- sqrt(b) * (log(time[i]) - mu - beta[1]*trt[i])
    S.norm[i] <- 1 - pnorm(epsilon.form.norm[i], 0, 1)
    f.norm[i] <- exp(-pow(epsilon.form.norm[i], 2) / 2) / sqrt(2*pi)
    h.norm[i] <- f.norm[i]/S.norm[i]  
  }
  
  ##priors for betas
  beta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1/2)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.logn.cus = file.path(dir, 'logn.cus.txt')
writeLines(modelstring.logn.cus, con=modelfile.logn.cus)

datlist.tau = list(N = nrow(df),
                   time = df$survT,
                   is.observed = df$dead,
                   trt = df$trt,
                   age = df$age,
                   sex = df$sex,
                   body = df$body,
                   tcic = df$tcic,
                   smoke = df$smoke,
                   pi = pi)

logn.cus = jags.parallel(data = datlist.tau, 
                         parameters.to.save = c('beta', 'mu', 'sigma'),
                         n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                         n.thin = 1, model.file = modelfile.logn.cus)
logn.cus$BUGSoutput[10]$summary[,c(1,3,7)]

freq.logn = flexsurvreg(formula = Surv(survT, dead) ~ I(trt), data = df, dist = 'lnorm')
freq.logn




modelstring.ll.cus = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
                       h.ll[i]*S.ll[i],
                       S.ll[i]))
    
    ##constructing hazard and Survival function
    
    beta_formula[i] <- beta[1] * trt[i] + 
                       beta[2] * age[i] + 
                       beta[3] * sex[i] +
                       beta[4] * body[i] +
                       beta[5] * tcic[i] + 
                       beta[6] * smoke[i] 
    
    epsilon.form.ll[i] <- b * (log(time[i]) - mu - beta_formula[i])
    S.ll[i] <- pow(1 + exp(epsilon.form.ll[i]), -1)
    h.ll[i] <- b/time[i] * pow(1 + exp(-epsilon.form.ll[i]), -1)
  }
  
  ##priors for betas
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  theta <- -mu/sigma
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.ll.cus = file.path(dir, 'll.cus.txt')
writeLines(modelstring.ll.cus, con=modelfile.ll.cus)

ll.cus = jags.parallel(data = datlist.tau, 
                       parameters.to.save = c('beta', 'theta', 'sigma'),
                       n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                       n.thin = 1, model.file = modelfile.ll.cus)
ll.cus$BUGSoutput[10]$summary[,c(1,3,7)]


####################################################################
## Fitting the models
####################################################################

modelstring.cusLN = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
                       h[i]*S[i],
                       S[i]))
    
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(alpha))) 
                                - mu - alpha - beta.form[i])
    S.tau[i] <- 1 - pnorm(epsilon.form.tau[i], 0, 1)
    f.tau[i] <- exp(-pow(epsilon.form.tau[i], 2) / 2) / sqrt(2*pi)
    h.tau[i] <- f.tau[i]/S.tau[i]
    
    ##otherwise
    epsilon.form.norm[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.norm[i] <- 1 - pnorm(epsilon.form.norm[i], 0, 1)
    f.norm[i] <- exp(-pow(epsilon.form.norm[i], 2) / 2) / sqrt(2*pi)
    h.norm[i] <- f.norm[i]/S.norm[i]  
    
    S[i] <- ifelse(trt[i],
                   S.norm[i],
                   ifelse(time[i] > tau, S.norm[i], S.tau[i]))
    h[i] <- ifelse(trt[i],
                   h.norm[i],
                   ifelse(time[i] > tau, h.norm[i], h.tau[i]))
  }
  
  ##priors for betas
  alpha ~ dnorm(0, 0.001)
  for(i in 1:5){
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0, 5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusLN = file.path(dir, 'cusLN.txt')
writeLines(modelstring.cusLN, con=modelfile.cusLN)

datlist.tau = list(N = nrow(df),
                   time = df$survT,
                   is.observed = df$dead,
                   trt = df$trt,
                   age = df$age,
                   sex = df$sex,
                   body = df$body,
                   tcic = df$tcic,
                   smoke = df$smoke,
                   pi = pi)

tau.ln.cus = jags.parallel(data = datlist.tau, 
                        parameters.to.save = c('beta', 'mu', 'sigma', 'tau'),
                        n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                        n.thin = 1, model.file = modelfile.cusLN)
##it does not like log(time[i] - tau - tau*exp(beta1)) to be log(negative number)

tau.ln.cus$BUGSoutput[10]$summary[,c(1,3,7)]

#traceplot(tau.cus)




modelstring.cusW = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(alpha))) 
                                - mu - alpha - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(alpha))) *
                exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:5){
    beta[i] ~ dnorm(0, 0.001)
  }
  alpha ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusW = file.path(dir, 'cusW.txt')
writeLines(modelstring.cusW, con=modelfile.cusW)

tau.wb.cus = jags.parallel(data = datlist.tau, 
                           parameters.to.save = c('alpha', 'beta', 'mu', 'sigma', 'tau'),
                           n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                           n.thin = 1, model.file = modelfile.cusW)
tau.wb.cus$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(tau.wb.cus)

set.seed(123)
auto.wb.tau = autorun.jags(model = modelfile.cusW, data = datlist.tau, 
                           monitor = c('alpha', 'beta', 'mu', 'sigma', 'tau'), 
                           n.chains = 3, startburnin = 5000, startsample = 6000, 
                           thin = 1, method = 'parallel')
auto.wb.tau


modelstring.cusLL = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(alpha))) 
                                - mu - alpha - beta.form[i])
    S.tau[i] <- pow(1 + exp(epsilon.form.tau[i]), -1)
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(alpha))) *
                pow(1 + exp(-epsilon.form.tau[i]), -1)
    
    ##otherwise
    epsilon.form.ll[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.ll[i] <- pow(1 + exp(epsilon.form.ll[i]), -1)
    h.ll[i] <- b/time[i] * pow(1 + exp(-epsilon.form.ll[i]), -1)
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.ll[i]),
                   S.ll[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.ll[i]),
                   h.ll[i])
  }
  
  ##priors for betas
  alpha ~ dnorm(0, 0.001)
  for(i in 1:5){
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  theta <- -mu/sigma
  tau ~ dunif(0.5, 5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusLL = file.path(dir, 'cusLL.txt')
writeLines(modelstring.cusLL, con=modelfile.cusLL)

tau.ll.cus = jags.parallel(data = datlist.tau, 
                           parameters.to.save = c('alpha', 'beta', 'mu', 'theta', 'sigma', 'tau'),
                           n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                           n.thin = 1, model.file = modelfile.cusLL)
tau.ll.cus$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(tau.ll.cus)

set.seed(123)
auto.ll.tau = autorun.jags(model = modelfile.cusLL, data = datlist.tau, 
                           monitor = c('alpha', 'beta', 'mu', 'sigma', 'tau'), 
                           n.chains = 3, startburnin = 5000, startsample = 6000, 
                           thin = 1, method = 'parallel')
auto.ll.tau




###########################################################################
###########################################################################
##changing the usual tau term into tau*exp(eta*Z)


modelstring.cusW.eta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
    
    eta.form[i] <- eta[1] * age[i] + 
                   eta[2] * sex[i] +
                   eta[3] * body[i] +
                   eta[4] * tcic[i] + 
                   eta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau + eta.form[i], S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau + eta.form[i], h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  beta[1] ~ dnorm(0.31950353, 100)
  beta[2] ~ dnorm(-0.03392898, 500)
  beta[3] ~ dnorm(0.17949740, 100)
  beta[4] ~ dnorm(-0.52902532, 100)
  beta[5] ~ dnorm(0.10020030, 100)
  beta[6] ~ dnorm(-0.07205756, 50)
  
  for(i in 1:5){
    eta[i] ~ dunif(-1.5, 1.5)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(2.94177717, 10)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(1.67411853, 2)I(0.5, )
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusW.eta = file.path(dir, 'cusW.eta.txt')
writeLines(modelstring.cusW.eta, con=modelfile.cusW.eta)

quantile(rnorm(842, mean = 0, sd = sqrt(1/1.25)),probs=c(.025,.975))

tau.wb.cus.eta = jags.parallel(data = datlist.tau, 
                               parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                               n.chains = 3, n.burnin = 8000, n.iter = 10000, 
                               n.thin = 1, model.file = modelfile.cusW.eta)
tau.wb.cus.eta$BUGSoutput[10]$summary[,c(1,3,7)]
gelman.diag(as.mcmc(tau.wb.cus.eta))
traceplot(tau.wb.cus.eta)

set.seed(123)
auto.wb.eta = autorun.jags(model = modelfile.cusW.eta, data = datlist.tau, 
                           monitor = c('beta', 'mu', 'sigma', 'tau', 'eta'), 
                           n.chains = 3, startburnin = 5000, startsample = 6000, 
                           thin = 1, method = 'parallel')
##93496 iterations, 1.1 hours taken


modelstring.cusLL.eta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
                    
    eta.form[i] <- eta[1] * age[i] + 
                   eta[2] * sex[i] +
                   eta[3] * body[i] +
                   eta[4] * tcic[i] + 
                   eta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)         
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau * exp(eta.form[i])
                              + tau * exp(beta[1] + eta.form[i])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- pow(1 + exp(epsilon.form.tau[i]), -1)
    h.tau[i] <- b / (abs(time[i] - tau * exp(eta.form[i]) 
                      + tau * exp(beta[1] + eta.form[i])))
                      * pow(1 + exp(-epsilon.form.tau[i]), -1)
    
    ##otherwise
    epsilon.form.ll[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.ll[i] <- pow(1 + exp(epsilon.form.ll[i]), -1)
    h.ll[i] <- b/time[i] * pow(1 + exp(-epsilon.form.ll[i]), -1)
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.ll[i]),
                   S.ll[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.ll[i]),
                   h.ll[i])
  }
  
  ##priors for betas
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  for(i in 1:5){
    eta[i] ~ dnorm(0, 10)
  }

  ##prior for mu and sigma of the log-logistic distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  theta <- -mu/sigma
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusLL.eta = file.path(dir, 'cusLL.eta.txt')
writeLines(modelstring.cusLL.eta, con=modelfile.cusLL.eta)

tau.ll.cus.eta = jags.parallel(data = datlist.tau, 
                               parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                               n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                               n.thin = 1, model.file = modelfile.cusLL.eta)
tau.ll.cus.eta$BUGSoutput[10]$summary[,c(1,3,7)]


time.cus = system.time(jags.parallel(data = datlist.tau, 
                   parameters.to.save = c('beta', 'tau', 'theta', 'sigma'),
                   n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                   n.thin = 1, model.file = modelfile.cusLL))
time.cus.tau = system.time(jags.parallel(data = datlist.tau, 
                   parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                   n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                   n.thin = 1, model.file = modelfile.cusLL.eta))
time.cus; time.cus.tau



###Changes:
##a. non-smoker set to 0, smoking or past smoker set to 1
##b. include a two-way interaction effect at T > tau between trt and each covariates. (with eta included)
##c. search for articles about Bayesian subgroup identification. (which subgroup has smaller tau, which subgroup has higher trt effect)

#https://journals.sagepub.com/doi/abs/10.1177/1740774510396933



modelstring.cusW.delta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(alpha
                              + delta.form[i]))) - mu - alpha - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(alpha + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  alpha ~ dnorm(0, 0.001)
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusW.delta = file.path(dir, 'cusW.delta.txt')
writeLines(modelstring.cusW.delta, con=modelfile.cusW.delta)

tau.wb.cus.delta = jags.parallel(data = datlist.tau, 
                      parameters.to.save = c('alpha', 'beta', 'mu', 'sigma', 'tau', 'delta'),
                      n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                      n.thin = 1, model.file = modelfile.cusW.delta)
tau.wb.cus.delta$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(tau.wb.cus.delta)

set.seed(123)
auto.wb.delta = autorun.jags(model = modelfile.cusW.delta, data = datlist.tau, 
                             monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                             n.chains = 3, startburnin = 5000, startsample = 6000, 
                             thin = 1, method = 'parallel')
auto.wb.delta


modelstring.cusLL.delta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau       
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(alpha + delta.form[i])))
                              - mu - alpha - delta.form[i] - beta.form[i])
    S.tau[i] <- pow(1 + exp(epsilon.form.tau[i]), -1)
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(alpha + delta.form[i])))
                    * pow(1 + exp(-epsilon.form.tau[i]), -1)
    
    ##otherwise
    epsilon.form.ll[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.ll[i] <- pow(1 + exp(epsilon.form.ll[i]), -1)
    h.ll[i] <- b/time[i] * pow(1 + exp(-epsilon.form.ll[i]), -1)
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.ll[i]),
                   S.ll[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.ll[i]),
                   h.ll[i])
  }
  
  ##priors for betas
  alpha ~ dnorm(0, 0.001)
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
    beta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the log-logistic distribution
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  theta <- -mu/sigma

  mu ~ dnorm(0, 0.001)
  tau ~ dunif(0.5, 5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
#dir = 'C:/Users/xingz/OneDrive/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusLL.delta = file.path(dir, 'cusLL.delta.txt')
writeLines(modelstring.cusLL.delta, con=modelfile.cusLL.delta)

tau.ll.cus.delta = 
  jags.parallel(data = datlist.tau, 
                parameters.to.save = c('alpha', 'beta', 'mu', 'sigma', 'tau', 'delta'),
                n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                n.thin = 1, model.file = modelfile.cusLL.delta)
tau.ll.cus.delta$BUGSoutput$summary[,c(1,3,7)]

set.seed(123)
auto.ll.delta = autorun.jags(model = modelfile.cusLL.delta, data = datlist.tau, 
                             monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                             n.chains = 3, startburnin = 5000, startsample = 6000, 
                             thin = 1, method = 'parallel')
auto.ll.delta

auto.ipAFT = autorun.jags(model = modelfile.cusLL.delta, data = datlist.tau, 
                          monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                          n.chains = 3, startburnin = 5000, startsample = 6000, 
                          thin = 1, method = 'parallel')
auto.ipAFT2 = autorun.jags(model = modelfile.cusLL.delta, data = datlist.tau, 
                           monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                           n.chains = 3, startburnin = 5000, startsample = 6000, 
                           thin = 1, method = 'parallel')


modelstring.cusW.delta.eta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    eta.form[i] <- eta[1] * age[i] + 
                   eta[2] * sex[i] +
                   eta[3] * body[i] +
                   eta[4] * tcic[i] + 
                   eta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau * exp(eta.form[i]) 
                              + tau * exp(beta[1] + delta.form[i] + eta.form[i])))
                              - mu - beta[1] - delta.form[i] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau * exp(eta.form[i])  
                    + tau * exp(beta[1] + delta.form[i] + eta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
    eta[i] ~ dnorm(0, 0.1)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusW.delta.eta = file.path(dir, 'cusW.delta.eta.txt')
writeLines(modelstring.cusW.delta.eta, con=modelfile.cusW.delta.eta)

tau.wb.cus.delta.eta = 
  jags.parallel(data = datlist.tau, 
    parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta', 'eta'),
    n.chains = 3, n.burnin = 6000, n.iter = 8000, 
    n.thin = 1, model.file = modelfile.cusW.delta.eta)
tau.wb.cus.delta.eta$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(tau.wb.cus.delta.eta)



###standardize age and include alpha and age only in the eta model.
###try excluding age




modelstring.cusLL.delta.eta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    eta.form[i] <- eta[1] * age[i] + 
                   eta[2] * sex[i] +
                   eta[3] * body[i] +
                   eta[4] * tcic[i] + 
                   eta[5] * smoke[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau * exp(eta.form[i]) 
                              + tau * exp(beta[1] + delta.form[i] + eta.form[i])))
                              - mu - beta[1] - delta.form[i] - beta.form[i])
    S.tau[i] <- pow(1 + exp(epsilon.form.tau[i]), -1)
    h.tau[i] <- b / (abs(time[i] - tau * exp(eta.form[i])  
                    + tau * exp(beta[1] + delta.form[i] + eta.form[i])))
                    * pow(1 + exp(-epsilon.form.tau[i]), -1)
    
    ##otherwise
    epsilon.form.ll[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.ll[i] <- pow(1 + exp(epsilon.form.ll[i]), -1)
    h.ll[i] <- b/time[i] * pow(1 + exp(-epsilon.form.ll[i]), -1)
    
    S[i] <- ifelse(trt[i],
                   S.ll[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.ll[i], S.tau[i]))
    h[i] <- ifelse(trt[i],
                   h.ll[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.ll[i], h.tau[i]))
  }
  
  ##priors for betas
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
    eta[i] ~ dnorm(0, 0.1)
  }

  ##prior for mu and sigma of the log-logistic distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  theta <- -mu/sigma
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusLL.delta.eta = file.path(dir, 'cusLL.delta.eta.txt')
writeLines(modelstring.cusLL.delta.eta, con=modelfile.cusLL.delta.eta)

tau.ll.cus.delta.eta = 
  jags.parallel(data = datlist.tau, 
                parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta', 'eta'),
                n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                n.thin = 1, model.file = modelfile.cusLL.delta.eta)
tau.ll.cus.delta.eta$BUGSoutput[10]$summary[,c(1,3,7)]
traceplot(tau.ll.cus.delta.eta)



##5 eta model with each covariate only.
##age, sex, body, tcic, smoke
modelstring.eta.baseline = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- 0
    
    eta.form[i] <- 0
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau * exp(eta.form[i]) 
                              + tau * exp(beta[1] + eta.form[i])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau * exp(eta.form[i]) 
                    + tau * exp(beta[1] + eta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  beta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.baseline = file.path(dir, 'eta.baseline.txt')
writeLines(modelstring.eta.baseline, con=modelfile.eta.baseline)

eta.baseline = jags.parallel(data = datlist.tau, 
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau'),
                             n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                             n.thin = 1, model.file = modelfile.eta.baseline)
eta.baseline$BUGSoutput[10]$summary[,c(1,3,7)]


modelstring.eta.age = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i]
    
    eta.form[i] <- eta[1] * age[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  eta[1] ~ dnorm(0, 10)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.age = file.path(dir, 'eta.age.txt')
writeLines(modelstring.eta.age, con=modelfile.eta.age)

eta.age = jags.parallel(data = datlist.tau, 
                        parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                        n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                        n.thin = 1, model.file = modelfile.eta.age)
eta.age$BUGSoutput[10]$summary[,c(1,3,7)]


modelstring.eta.sex = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * sex[i]
    
    eta.form[i] <- eta[1] * sex[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  eta[1] ~ dnorm(0, 10)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.sex = file.path(dir, 'eta.sex.txt')
writeLines(modelstring.eta.sex, con=modelfile.eta.sex)

eta.sex = jags.parallel(data = datlist.tau, 
                        parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                        n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                        n.thin = 1, model.file = modelfile.eta.sex)
eta.sex$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.eta.body = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * body[i]
    
    eta.form[i] <- eta[1] * body[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  eta[1] ~ dnorm(0, 10)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.body = file.path(dir, 'eta.body.txt')
writeLines(modelstring.eta.body, con=modelfile.eta.body)

eta.body = jags.parallel(data = datlist.tau, 
                        parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                        n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                        n.thin = 1, model.file = modelfile.eta.body)
eta.body$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.eta.tcic = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * tcic[i]
    
    eta.form[i] <- eta[1] * tcic[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  eta[1] ~ dnorm(0, 10)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.tcic = file.path(dir, 'eta.tcic.txt')
writeLines(modelstring.eta.tcic, con=modelfile.eta.tcic)

eta.tcic = jags.parallel(data = datlist.tau, 
                        parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                        n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                        n.thin = 1, model.file = modelfile.eta.tcic)
eta.tcic$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.eta.smoke = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * smoke[i]
    
    eta.form[i] <- eta[1] * smoke[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                              - mu - beta[1] - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(beta[1]) + eta.form[i] * exp(beta[1])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta[1]*trt[i] - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau*exp(eta.form[i]), h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  eta[1] ~ dnorm(0, 10)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.smoke = file.path(dir, 'eta.smoke.txt')
writeLines(modelstring.eta.smoke, con=modelfile.eta.smoke)

eta.smoke = jags.parallel(data = datlist.tau, 
                        parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta'),
                        n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                        n.thin = 1, model.file = modelfile.eta.smoke)
eta.smoke$BUGSoutput[10]$summary[,c(1,3,7)]

#tau.wb.cus$BUGSoutput[10]$summary[,c(1,3,7)]
#tau.wb.cus.eta$BUGSoutput[10]$summary[,c(1,3,7)]
cbind(tau.wb.cus$BUGSoutput[10]$summary[,1], tau.wb.cus.eta$BUGSoutput[10]$summary[,1])

eta.baseline$BUGSoutput[10]$summary[,c(1,3,7)]
eta.age$BUGSoutput[10]$summary[,c(1,3,7)]
eta.sex$BUGSoutput[10]$summary[,c(1,3,7)]
eta.body$BUGSoutput[10]$summary[,c(1,3,7)]
eta.tcic$BUGSoutput[10]$summary[,c(1,3,7)]
eta.smoke$BUGSoutput[10]$summary[,c(1,3,7)]

eta.smoke.limit$BUGSoutput[10]$summary[,c(1,3,7)]

traceplot(eta.baseline)
traceplot(eta.age) ##tau's prior is not good
traceplot(eta.sex)##beta1, eta, tau
traceplot(eta.body)##beta1, beta2, eta, mu, sigma, tau
traceplot(eta.tcic)##beta1, eta, mu, tau
traceplot(eta.smoke)##beta1, eta, mu, sigma, tau

##traceplot indicates that models with eta has multimodality.

xyplot(as.mcmc(eta.age))
densityplot(as.mcmc(eta.baseline))




########################################################################
########################################################################


modelstring.cusW.delta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:6){
    beta[i] ~ dnorm(0, 0.001)
  }
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusW.delta = file.path(dir, 'cusW.delta.txt')
writeLines(modelstring.cusW.delta, con=modelfile.cusW.delta)

tau.wb.cus.delta = jags.parallel(data = datlist.tau, 
                                 parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                                 n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                                 n.thin = 1, model.file = modelfile.cusW.delta)
tau.wb.cus.delta$BUGSoutput[10]$summary[,c(1,3,7)]
traceplot(tau.wb.cus.delta)

################################################################################

modelstring.delta.age = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] 
    
    delta.form[i] <- delta[1] * age[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.age = file.path(dir, 'delta.age.txt')
writeLines(modelstring.delta.age, con=modelfile.delta.age)

delta.age = jags.parallel(data = datlist.tau, 
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                          n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                          n.thin = 1, model.file = modelfile.delta.age)
delta.age$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.delta.sex = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * sex[i] 
    
    delta.form[i] <- delta[1] * sex[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.sex = file.path(dir, 'delta.sex.txt')
writeLines(modelstring.delta.sex, con=modelfile.delta.sex)

delta.sex = jags.parallel(data = datlist.tau, 
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                          n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                          n.thin = 1, model.file = modelfile.delta.sex)
delta.sex$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.delta.body = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * body[i] 
    
    delta.form[i] <- delta[1] * body[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.body = file.path(dir, 'delta.body.txt')
writeLines(modelstring.delta.body, con=modelfile.delta.body)

delta.body = jags.parallel(data = datlist.tau, 
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                          n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                          n.thin = 1, model.file = modelfile.delta.body)
delta.body$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.delta.tcic = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * tcic[i] 
    
    delta.form[i] <- delta[1] * tcic[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.tcic = file.path(dir, 'delta.tcic.txt')
writeLines(modelstring.delta.tcic, con=modelfile.delta.tcic)

delta.tcic = jags.parallel(data = datlist.tau, 
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                          n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                          n.thin = 1, model.file = modelfile.delta.tcic)
delta.tcic$BUGSoutput[10]$summary[,c(1,3,7)]



modelstring.delta.smoke = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * smoke[i] 
    
    delta.form[i] <- delta[1] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.smoke = file.path(dir, 'delta.smoke.txt')
writeLines(modelstring.delta.smoke, con=modelfile.delta.smoke)

delta.smoke = jags.parallel(data = datlist.tau, 
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                          n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                          n.thin = 1, model.file = modelfile.delta.smoke)
delta.smoke$BUGSoutput[10]$summary[,c(1,3,7)]



delta.age$BUGSoutput[10]$summary[,c(1,3,7)]
delta.sex$BUGSoutput[10]$summary[,c(1,3,7)]
delta.body$BUGSoutput[10]$summary[,c(1,3,7)]
delta.tcic$BUGSoutput[10]$summary[,c(1,3,7)]
delta.smoke$BUGSoutput[10]$summary[,c(1,3,7)]

#traceplot(delta.smoke)





modelstring.delta.smoke2 = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * smoke[i] 
    
    delta.form[i] <- delta[1] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i]-delta.form[i]-beta[1]*trt[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.smoke2 = file.path(dir, 'delta.smoke2.txt')
writeLines(modelstring.delta.smoke2, con=modelfile.delta.smoke2)

delta.smoke2 = jags.parallel(data = datlist.tau, 
                            parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                            n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                            n.thin = 1, model.file = modelfile.delta.smoke2)
delta.smoke2$BUGSoutput[10]$summary[,c(1,3,7)]


modelstring.delta.body2 = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * body[i] 
    
    delta.form[i] <- delta[1] * body[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i]-delta.form[i]-beta[1]*trt[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  
    delta[1] ~ dnorm(0, 0.001)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.body2 = file.path(dir, 'delta.body2.txt')
writeLines(modelstring.delta.body2, con=modelfile.delta.body2)

delta.body2 = jags.parallel(data = datlist.tau, 
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                             n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                             n.thin = 1, model.file = modelfile.delta.body2)
delta.body2$BUGSoutput[10]$summary[,c(1,3,7)]
traceplot(delta.body2)




#P(time < tau | covariates)
#P(time < tau + eta*Z | covariates)

#get result from simple tau model, and then limit the prior for the delta model
tau.res = as.data.frame(tau.wb.cus$BUGSoutput[10]$summary[,c(1,2)])
tau.res$prec = 1/(tau.res$sd)^2
tau.res


modelstring.delta.mod = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  beta[1] ~ dnorm(0.31950353, 157.3212628)
  beta[2] ~ dnorm(-0.03392898, 787.5179331)
  beta[3] ~ dnorm(0.17949740, 155.4812210)
  beta[4] ~ dnorm(-0.52902532, 156.6206606)
  beta[5] ~ dnorm(0.10020030, 188.3348114)
  beta[6] ~ dnorm(-0.07205756, 99.8025067)
  
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(2.94177717, 66.3486655)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.mod = file.path(dir, 'delta.mod.txt')
writeLines(modelstring.delta.mod, con=modelfile.delta.mod)

delta.mod = jags.parallel(data = datlist.tau, 
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                          n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                          n.thin = 1, model.file = modelfile.delta.mod)
delta.mod$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(delta.mod)


modelstring.delta.mod2 = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[2] * age[i] + 
                    beta[3] * sex[i] +
                    beta[4] * body[i] +
                    beta[5] * tcic[i] + 
                    beta[6] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau*exp(eta*Z)        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(beta[1]
                              + delta.form[i]))) - mu - beta[1] - delta.form[i]
                              - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(beta[1] + delta.form[i])))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  beta[1] ~ dnorm(0.31950353, 50)
  beta[2] ~ dnorm(-0.03392898, 100)
  beta[3] ~ dnorm(0.17949740, 50)
  beta[4] ~ dnorm(-0.52902532, 50)
  beta[5] ~ dnorm(0.10020030, 50)
  beta[6] ~ dnorm(-0.07205756, 30)
  
  
  for(i in 1:5){
    delta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(2.94177717, 20)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5, 3)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.mod2 = file.path(dir, 'delta.mod2.txt')
writeLines(modelstring.delta.mod2, con=modelfile.delta.mod2)

delta.mod2 = jags.parallel(data = datlist.tau, 
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta'),
                           n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                           n.thin = 1, model.file = modelfile.delta.mod2)
delta.mod2$BUGSoutput[10]$summary[,c(1,3,7)]
traceplot(delta.mod)


########################################################################
########################################################################
########################################################################
########################################################################
X = df[, -c(1, 2, 9:12)]

##calculating the mean of beta^T X from all mcmc iterations
tau.chain = rbind(auto.ll.delta$mcmc[[1]], auto.ll.delta$mcmc[[2]], auto.ll.delta$mcmc[[3]])
beta.X = tau.chain[,1:6] %*% t(X)
dim(beta.X)

prob.temp = matrix(nrow = nrow(beta.X), ncol = ncol(beta.X))
for(i in 1:ncol(beta.X)){
  prob.temp[, i] = 1 - 1/(1 + exp(
    1/tau.chain[,13]*(log(tau.chain[,14])-tau.chain[,12]-beta.X[,i])))
}
##probability of patients' survival time being more than tau for 850 patients.
prob.morethantau = apply(1 - prob.temp, 2, FUN = mean)




########################################################################

##calculate the treatment effect of each patient from the delta model
##exp(alpha+delta*X)
X.trt1 = X
X.trt1$trt = 1

delta.chain = rbind(auto.ll.delta$mcmc[[1]], 
                    auto.ll.delta$mcmc[[2]], auto.ll.delta$mcmc[[3]])
delta.X = delta.chain[,c(1, 7:11)] %*% t(as.matrix(X.trt1))
prob.temp2 = matrix(data = 0, nrow = nrow(delta.X), ncol = ncol(delta.X))
for(i in 1:ncol(delta.X)){
    prob.temp2[,i] = exp(delta.X[,i])
}

delta.treat = apply(prob.temp2, 2, FUN = mean)


delta.chain = rbind(auto.ll.delta$mcmc[[1]], 
                    auto.ll.delta$mcmc[[2]], auto.ll.delta$mcmc[[3]])
delta.X = delta.chain[,c(1, 7:11)] %*% t(as.matrix(X.trt1))
prob.temp2 = matrix(data = 0, nrow = nrow(delta.X), ncol = ncol(delta.X))
delta.treat = apply(delta.X, 2, FUN = mean)
delta.treat = exp(delta.treat)
########################################################################

X[,-2] = apply(X[,-2], 2, function(x) as.factor(x))
Xplot = cbind(X, prob.morethantau, delta.treat)
Xplot$elderly = as.factor(ifelse(X$age*sd(data$BAGE)+mean(data$BAGE)>=64.5, 1, 0))

library(ggpubr)

plot.x1 = ggplot(Xplot,aes(x=prob.morethantau, color=elderly)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'PRT', color = 'Age > 65')
plot.x2 = ggplot(Xplot,aes(x=prob.morethantau, color=sex)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'PRT', color = 'Sex')
plot.x3 = ggplot(Xplot,aes(x=prob.morethantau, color=body)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'PRT', color = 'Body')
plot.x4 = ggplot(Xplot,aes(x=prob.morethantau, color=tcic)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'PRT', color = 'TCIC')
plot.x5 = ggplot(Xplot,aes(x=prob.morethantau, color=smoke)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'PRT', color = 'Smoke')
plot.x = ggarrange(plot.x1, plot.x2, plot.x3, plot.x4, plot.x5, ncol=3, nrow=2)
annotate_figure(plot.x)



plot.y1 = ggplot(Xplot, aes(x=delta.treat, color=elderly)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'TE', color = 'Age > 65')
plot.y2 = ggplot(Xplot, aes(x=delta.treat, color=sex)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'TE', color = 'Sex')
plot.y3 = ggplot(Xplot, aes(x=delta.treat, color=body)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'TE', color = 'Body')
plot.y4 = ggplot(Xplot, aes(x=delta.treat, color=tcic)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'TE', color = 'TCIC')
plot.y5 = ggplot(Xplot, aes(x=delta.treat, color=smoke)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill='white', position = 'dodge') + 
  geom_density() + 
  labs(x = 'TE', color = 'Smoke')

plot.y = ggarrange(plot.y1, plot.y2, plot.y3, plot.y4, plot.y5, ncol=3, nrow=2)
annotate_figure(plot.y)

library(viridis)
##by age
plot.xy1 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=age)) + 
  scale_color_viridis(option = "plasma", direction = -1) +
  labs(x = 'PRT', y = 'TE', color = 'Age') 
##by sex
plot.xy2 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=sex)) + 
  labs(x = 'PRT', y = 'TE', color = 'Sex') 
##by body
plot.xy3 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=body)) + 
  labs(x = 'PRT', y = 'TE', color = 'Body')
##by tcic
plot.xy4 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=tcic)) + 
  labs(x = 'PRT', y = 'TE', color = 'TCIC')
##by smoke
plot.xy5 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=smoke)) + 
  labs(x = 'PRT', y = 'TE', color = 'Smoke')
plot.xy = ggarrange(plot.xy1, plot.xy2, plot.xy3, plot.xy4, plot.xy5, ncol=3, nrow=2)
annotate_figure(plot.xy)

ggsave(
  "trt vs prob.png",
  width = 10.5,
  height = 5,
  dpi = 800
)

#############################################################################
#############################################################################

#wb.logl = autorun.jags(model = modelfile.weibull.cus, data = datlist.tau, 
#                       monitor = c('l'), n.chains = 3, startburnin = 5000, 
#                       startsample = 6000, thin = 1, method = 'parallel')
#ll.logl = autorun.jags(model = modelfile.ll.cus, data = datlist.tau, 
#                       monitor = c('l'), n.chains = 3, startburnin = 5000, 
#                       startsample = 6000, thin = 1, method = 'parallel')

#tau.wb.logl = autorun.jags(model = modelfile.cusW, data = datlist.tau, 
#                        monitor = c('l'), n.chains = 3, startburnin = 5000, 
#                        startsample = 6000, thin = 1, method = 'parallel')
#delta.wb.logl = autorun.jags(model = modelfile.cusW.delta, data = datlist.tau, 
#                          monitor = c('l'), n.chains = 3, startburnin = 5000, 
#                          startsample = 6000, thin = 1, method = 'parallel')
#tau.ll.logl = autorun.jags(model = modelfile.cusLL, data = datlist.tau, 
#                           monitor = c('l'), n.chains = 3, startburnin = 5000, 
#                           startsample = 6000, thin = 1, method = 'parallel')
#delta.ll.logl = autorun.jags(model = modelfile.cusLL.delta, data = datlist.tau, 
#                             monitor = c('l'), n.chains = 3, startburnin = 5000, 
#                              startsample = 6000, thin = 1, method = 'parallel')



wb.logl.chain = rbind(wb.logl$mcmc[[1]], wb.logl$mcmc[[2]], 
                          wb.logl$mcmc[[3]])
wb.logl.mean = apply(wb.logl.chain, 2, function(x) log(mean(exp(x))))
wb.lppd = sum(wb.logl.mean)
wb.logl.var = apply(wb.logl.chain, 2, var)
wb.p_waic = sum(wb.logl.var)
(wb.waic = -2 * (wb.lppd - wb.p_waic))

ll.logl.chain = rbind(ll.logl$mcmc[[1]], ll.logl$mcmc[[2]], 
                      ll.logl$mcmc[[3]])
ll.logl.mean = apply(ll.logl.chain, 2, function(x) log(mean(exp(x))))
ll.lppd = sum(ll.logl.mean)
ll.logl.var = apply(ll.logl.chain, 2, var)
ll.p_waic = sum(ll.logl.var)
(ll.waic = -2 * (ll.lppd - ll.p_waic))

tau.wb.logl.chain = rbind(tau.wb.logl$mcmc[[1]], tau.wb.logl$mcmc[[2]], 
                          tau.wb.logl$mcmc[[3]])
tau.wb.logl.mean = apply(tau.wb.logl.chain, 2, function(x) log(mean(exp(x))))
tau.wb.lppd = sum(tau.wb.logl.mean)
tau.wb.logl.var = apply(tau.wb.logl.chain, 2, var)
tau.wb.p_waic = sum(tau.wb.logl.var)
(tau.wb.waic = -2 * (tau.wb.lppd - tau.wb.p_waic))

tau.ll.logl.chain = rbind(tau.ll.logl$mcmc[[1]], tau.ll.logl$mcmc[[2]], 
                          tau.ll.logl$mcmc[[3]])
tau.ll.logl.mean = apply(tau.ll.logl.chain, 2, function(x) log(mean(exp(x))))
tau.ll.lppd = sum(tau.ll.logl.mean)
tau.ll.logl.var = apply(tau.ll.logl.chain, 2, var)
tau.ll.p_waic = sum(tau.ll.logl.var)
(tau.ll.waic = -2 * (tau.ll.lppd - tau.ll.p_waic))

delta.wb.logl.chain = rbind(delta.wb.logl$mcmc[[1]], delta.wb.logl$mcmc[[2]], 
                         delta.wb.logl$mcmc[[3]])
delta.wb.logl.mean = apply(delta.wb.logl.chain, 2, function(x) log(mean(exp(x))))
delta.wb.lppd = sum(delta.wb.logl.mean)
delta.wb.logl.var = apply(delta.wb.logl.chain, 2, var)
delta.wb.p_waic = sum(delta.wb.logl.var)
(delta.wb.waic = -2 * (delta.wb.lppd - delta.wb.p_waic))

delta.ll.logl.chain = rbind(delta.ll.logl$mcmc[[1]], delta.ll.logl$mcmc[[2]], 
                            delta.ll.logl$mcmc[[3]])
delta.ll.logl.mean = apply(delta.ll.logl.chain, 2, function(x) log(mean(exp(x))))
delta.ll.lppd = sum(delta.ll.logl.mean)
delta.ll.logl.var = apply(delta.ll.logl.chain, 2, var)
delta.ll.p_waic = sum(delta.ll.logl.var)
(delta.ll.waic = -2 * (delta.ll.lppd - delta.ll.p_waic))


waic = c(wb.waic, ll.waic, tau.wb.waic, tau.ll.waic, delta.wb.waic, delta.ll.waic)
names(waic) = c('wb', 'll', 'wb.tau', 'll.tau', 'wb.delta', 'll.delta')
waic

waic.order = c(delta.ll.waic, tau.ll.waic, ll.waic, delta.wb.waic, wb.waic, tau.wb.waic)
names(waic.order) = c('ll.delta', 'll.tau', 'll', 'wb.delta', 'wb', 'wb.tau')
waic.order





##forest plot

library(forestplot)

ci.morethantau = cbind(apply(1- prob.temp, 2, FUN = function(x){quantile(x, probs=c(0.025))}),
                       apply(1- prob.temp, 2, FUN = function(x){quantile(x, probs=c(0.975))}))
##get the quantile of 2.5% and 97.5%
ci.name = as.character(as.vector(1:842))
ci.plot = as.data.frame(cbind(ci.name, prob.morethantau, ci.morethantau))
library(dplyr)
ci.plot <- ci.plot %>%
  arrange(prob.morethantau)
forestplot(mean=ci.plot[,2], lower=ci.plot[,3], 
           upper=ci.plot[,4], labeltext=ci.plot[,1], clip=c(0.725,Inf),
           boxsize=0.02, xlim = c(0.725, 1), xticks = seq(0.725, 1, by = 0.025),
           zero = 0.725, lineheight = "auto")



library(Hmisc)
find.q = function(x){
  return(cbind(mean(x),
               quantile(x,probs=(0.025)),
               quantile(x,probs=(0.975))))
}

prob.plot = rbind(find.q(prob.morethantau[Xplot$elderly==0]),
                  find.q(prob.morethantau[Xplot$elderly==1]),
                  find.q(prob.morethantau[X$sex==0]),
                  find.q(prob.morethantau[X$sex==1]),
                  find.q(prob.morethantau[X$smoke==0]),
                  find.q(prob.morethantau[X$smoke==1]),
                  find.q(prob.morethantau[X$body==0]),
                  find.q(prob.morethantau[X$body==1]),
                  find.q(prob.morethantau[X$tcic==0]),
                  find.q(prob.morethantau[X$tcic==1]))
rownames(prob.plot) = c('Age<65', 'Age>65', 'male', 'female', 'no smoke', 'smoke',
                     'body=0', 'body=1', 'tcic=0', 'tcic=1')
prob.plot.data = tibble::tibble(mean = prob.plot[,1],
                           lower = prob.plot[,2],
                           upper = prob.plot[,3],
                           factor = rownames(prob.plot))
prob.plot.data |>
  forestplot(labeltext = c(factor),
             zero = 0.85,
             xlog = FALSE,
             boxsize = 0.1,          # Increase box size
             lwd.ci = 2,             # Thicker confidence interval lines
             lwd.zero = 2,           # Thicker zero line
             cex = 1.2,              # Increase text size
             xlim = c(0.85, 1),       # Set the x-axis limits to focus on the relevant range
             xticks = seq(0.85, 1, by = 0.025),  # Adjust the tick marks
             col = forestplot::fpColors(box = "red", line = "darkred", summary = "red"),  # Use high-contrast colors
             )|>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue")|>
  fp_set_zebra_style("#EFEFEF")



delta.plot = rbind(find.q(delta.treat[Xplot$elderly==0]),
                   find.q(delta.treat[Xplot$elderly==1]),
                   find.q(delta.treat[X$sex==0]),
                   find.q(delta.treat[X$sex==1]),
                   find.q(delta.treat[X$smoke==0]),
                   find.q(delta.treat[X$smoke==1]),
                   find.q(delta.treat[X$body==0]),
                   find.q(delta.treat[X$body==1]),
                   find.q(delta.treat[X$tcic==0]),
                   find.q(delta.treat[X$tcic==1]))
rownames(delta.plot) = c('Age<65', 'Age>65', 'male', 'female', 'no smoke', 'smoke',
                        'body=0', 'body=1', 'tcic=0', 'tcic=1')
delta.plot.data = tibble::tibble(mean = delta.plot[,1],
                           lower = delta.plot[,2],
                           upper = delta.plot[,3],
                           factor = rownames(delta.plot))
delta.plot.data |>
  forestplot(labeltext = c(factor),
             zero = 0.9,
             xlog = FALSE,
             boxsize = 0.1,          # Increase box size
             lwd.ci = 2,             # Thicker confidence interval lines
             lwd.zero = 2,           # Thicker zero line
             cex = 1.2,              # Increase text size
             xlim = c(0.9, 2.6),       # Set the x-axis limits to focus on the relevant range
             xticks = seq(0.9, 2.6, by = 0.1),  # Adjust the tick marks
             col = forestplot::fpColors(box = "red", line = "darkred", summary = "red"),  # Use high-contrast colors
             title = "Treatment Effect")|>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue")|>
  fp_set_zebra_style("#EFEFEF")

load('.RData')
#(x-min(x))/(max(x)-min(x))
library(dendextend)
library(pheatmap)
cluster.plot = cbind(prob.morethantau, delta.treat)
colnames(cluster.plot) = c('PRT', 'TE')
cluster.std = apply(cluster.plot, 2, FUN=function(x){return((x-mean(x))/sd(x))})
dist.matrix = dist(cluster.std)
hc <- hclust(dist.matrix, method = 'average')
clusters = cutree(hc, k = 5)
plot(cluster.plot, col = clusters, pch = 19, cex=0.5, xlab = "PRT", 
     ylab = "TE")

plot(hc, main = "Cluster Dendrogram", xlab = "Samples", sub = "", cex = 0.9)
rect.hclust(hc, k = 3, border = c("red", "green", "blue", "cyan"))

##Elbow method 
wss <- function(k) {
  clusters <- cutree(hc, k)
  sum(sapply(unique(clusters), function(i) sum(dist.matrix[clusters == i])^2))
}
k_values <- 2:10  # Test different cluster values
wss_values <- sapply(k_values, wss)
plot(k_values, wss_values, type = "b", pch = 19,
     xlab = "Number of Clusters (k)", ylab = "Within Sum of Squares",
     main = "Elbow Method for Optimal Clusters")

##Silhouette method
silhouette_score <- function(k) {
  clusters <- cutree(hc, k)
  ss <- silhouette(clusters, dist.matrix)
  return(mean(ss[, 3]))
}
sil_scores <- sapply(k_values, silhouette_score)
plot(k_values, sil_scores, type = "b", pch = 19,
     xlab = "Number of Clusters (k)", ylab = "Average Silhouette Score",
     main = "Silhouette Method for Optimal Clusters")


cluster_colors <- c("#F8766D",  # Red-like
                    "#00BFC4",  # Blue-like
                    "#7CAE00",  # Dark Green-like
                    "#C77CFF",  # Purple-like
                    "#FFA500")  # Orange-like

Cluster2 = as.factor(clusters)
cluster_analysis2 <- ggplot(as.data.frame(cluster.plot), aes(x = `PRT`, 
                                        y = `TE`, 
                                        color = Cluster2)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cluster_colors[2:1]) +
  labs(x = "PRT", y = "TE", color = '2 Clusters')

##3 clusters
clusters3 = cutree(hc, k = 3)
Cluster3 = as.factor(clusters3)
cluster_analysis3 <- ggplot(as.data.frame(cluster.plot), aes(x = `PRT`, 
                                                            y = `TE`, 
                                                            color = Cluster3)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cluster_colors[c(2,4,1)]) +
  labs(x = "PRT", y = "TE", color = '3 Clusters')

##4 clusters
clusters4 = cutree(hc, k = 4)
Cluster4 = as.factor(clusters4)
cluster_analysis4 <- ggplot(as.data.frame(cluster.plot), aes(x = `PRT`, 
                                                             y = `TE`, 
                                                             color = Cluster4)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cluster_colors[c(2,4,3,1)]) +
  labs(x = "PRT", y = "TE", color = '4 Clusters')

##5 clusters
clusters5 = cutree(hc, k = 5)
##rearrange the order of the cluster to make the plot more easier to examine
mapping <- c(`3` = 1, `4` = 2, `5` = 3, `2` = 4, `1` = 5)
clusters5 = unname(mapping[as.character(clusters5)])
Cluster5 = as.factor(clusters5)
cluster_analysis5 <- ggplot(as.data.frame(cluster.plot), aes(x = `PRT`, 
                                                             y = `TE`, 
                                                             color = Cluster5)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cluster_colors[c(3,5,1,4,2)]) +
  labs(x = "PRT", y = "TE", color = '5 Clusters')

plot_cluster = ggarrange(cluster_analysis2, cluster_analysis3, 
                         cluster_analysis4, cluster_analysis5, ncol=2, nrow=2)
annotate_figure(plot_cluster)

plot.xy = ggarrange(plot.xy1, plot.xy2, plot.xy3, plot.xy4, plot.xy5, cluster_analysis5, 
                    ncol=3, nrow=2)
annotate_figure(plot.xy)

pheatmap(t(cluster.std),
         cluster_cols = hc,
         cluster_rows = F,
         clustering_method = 'ward.D2',
         cutree_cols = 5,
         treeheight_col = 100,
         annotation_names_row = 1,
         annotation_colors = 'red')
dev.new()
pheatmap(t(cluster.std),
         color = colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))(100),
         cluster_cols = hc,
         cluster_rows = F,
         clustering_method = 'average',
         cutree_cols = 5,
         treeheight_col = 100,
         annotation_names_row = 1
)
png('KM oak.png', width=7, height=5, units='in', res=800)
dev.off()




##check the characteristics of each clusters
c1 =lapply(Xplot[clusters5==1,], as.numeric)
c2 =lapply(Xplot[clusters5==2,], as.numeric)
c3 =lapply(Xplot[clusters5==3,], as.numeric)
c4 =lapply(Xplot[clusters5==4,], as.numeric)
c5 =lapply(Xplot[clusters5==5,], as.numeric)

##cluster characteristics
cluster.char.temp = matrix(
        cbind(lapply(c1, FUN=sum), lapply(c2, FUN=sum), 
                     lapply(c3, FUN=sum), lapply(c4, FUN=sum), 
                     lapply(c5, FUN=sum)), nrow=9, ncol=5)
cluster.char.temp <- matrix(as.numeric(cluster.char.temp), nrow = nrow(cluster.char.temp),
                            ncol = ncol(cluster.char.temp))
cluster.char = matrix()
c.length = c(length(c1$sex),length(c2$sex),length(c3$sex),length(c4$sex),length(c5$sex))
##age
cluster.char = rbind(c(mean(c1$age)*sd(data$BAGE)+mean(data$BAGE),
                       mean(c2$age)*sd(data$BAGE)+mean(data$BAGE),
                       mean(c3$age)*sd(data$BAGE)+mean(data$BAGE),
                       mean(c4$age)*sd(data$BAGE)+mean(data$BAGE),
                       mean(c5$age)*sd(data$BAGE)+mean(data$BAGE)))
##male and female
cluster.char = rbind(cluster.char, c.length-cluster.char.temp[3,])
cluster.char = rbind(cluster.char, cluster.char.temp[3,])
##body
cluster.char = rbind(cluster.char, c.length-cluster.char.temp[5,])
cluster.char = rbind(cluster.char, cluster.char.temp[5,])
##tcic
cluster.char = rbind(cluster.char, c.length-cluster.char.temp[6,])
cluster.char = rbind(cluster.char, cluster.char.temp[6,])
##smoke
cluster.char = rbind(cluster.char, c.length-cluster.char.temp[4,])
cluster.char = rbind(cluster.char, cluster.char.temp[4,])
##PRT and TE
cluster.char = rbind(cluster.char, c(mean(c1$prob.morethantau),
                                     mean(c2$prob.morethantau),
                                     mean(c3$prob.morethantau),
                                     mean(c4$prob.morethantau),
                                     mean(c5$prob.morethantau)))
cluster.char = rbind(cluster.char, c(mean(c1$delta.treat),
                                     mean(c2$delta.treat),
                                     mean(c3$delta.treat),
                                     mean(c4$delta.treat),
                                     mean(c5$delta.treat)))

colnames(cluster.char) = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4',
                           'Cluster 5')
rownames(cluster.char) = c('Age', 'Male', 'Female', 'body=0', 'body=1','tcic=0',
                           'tcic=1', 'non-smoker', 'smoker', 'PRT', 'TE')
cluster.char


##quantile of aPRT and aTE
q.prt = rbind(quantile(c1$prob.morethantau, c(0.25, 0.75)),
              quantile(c2$prob.morethantau, c(0.25, 0.75)),
              quantile(c3$prob.morethantau, c(0.25, 0.75)),
              quantile(c4$prob.morethantau, c(0.25, 0.75)),
              quantile(c5$prob.morethantau, c(0.25, 0.75)))
q.te = rbind(quantile(c1$delta.treat, c(0.25, 0.75)),
             quantile(c2$delta.treat, c(0.25, 0.75)),
             quantile(c3$delta.treat, c(0.25, 0.75)),
             quantile(c4$delta.treat, c(0.25, 0.75)),
             quantile(c5$delta.treat, c(0.25, 0.75)))



##averaged treatment effect of ipAFT model
exp(mean(log(delta.treat)))

##characteristics of patients with delta.treat <= 1
lapply(lapply(Xplot[delta.treat<=1,], as.numeric), FUN=mean)







modelstring.cusW.eta = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    eta.form[i] <- tau^exp( 
                   eta[1] * age[i] + 
                   eta[2] * sex[i] +
                   eta[3] * body[i] +
                   eta[4] * tcic[i] + 
                   eta[5] * smoke[i])
    
    ##when treatment = T and Time_i >= tau + eta*W        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - eta.form[i]
                              + exp(alpha)*eta.form[i]))
                              - mu - alpha - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - eta.form[i]
                    + exp(alpha)*eta.form[i]))
                    * exp(epsilon.form.tau[i])
    
    ##otherwise
    epsilon.form.wb[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.wb[i] <- exp(-exp(epsilon.form.wb[i]))
    h.wb[i] <- (b/time[i]) * exp(epsilon.form.wb[i])
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > eta.form[i], S.tau[i], S.wb[i]),
                   S.wb[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > eta.form[i], h.tau[i], h.wb[i]),
                   h.wb[i])
  }
  
  ##priors for betas
  #alpha ~ dnorm(0, 0.01)
  #for(i in 1:5){
  #  beta[i] ~ dnorm(0, 0.01)
  #}
  
  alpha ~ dnorm(0.31950353, 100)
  beta[1] ~ dnorm(-0.03392898, 500)
  beta[2] ~ dnorm(0.17949740, 100)
  beta[3] ~ dnorm(-0.52902532, 100)
  beta[4] ~ dnorm(0.10020030, 100)
  beta[5] ~ dnorm(-0.07205756, 50)
  
  for(i in 1:5){
    eta[i] ~ dnorm(0, 1/2)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.01)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(1.6835, 5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
#dir = 'C:/Users/xingz/OneDrive/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusW.eta = file.path(dir, 'cusW.eta.txt')
writeLines(modelstring.cusW.eta, con=modelfile.cusW.eta)

tau.wb.cus.eta = jags.parallel(data = datlist.tau, 
                               parameters.to.save = c('alpha', 'beta', 'mu', 'sigma', 'tau', 'eta'),
                               n.chains = 3, n.burnin = 8000, n.iter = 10000, 
                               n.thin = 1, model.file = modelfile.cusW.eta)
tau.wb.cus.eta$BUGSoutput[10]$summary[,c(1,3,7)]

auto.eta = autorun.jags(model = modelfile.cusW.eta, data = datlist.tau, 
                        monitor = c('alpha', 'beta', 'mu', 'sigma', 'tau', 'eta'), 
                        n.chains = 3, startburnin = 5000, startsample = 6000, 
                        thin = 1, method = 'parallel')
auto.eta




modelstring.cusLL.delta2 = '
data{
  for(z in 1:N){
    zeros[z] <- 0
  }
  C <- 10000
}

model{
  for(i in 1 : N) {
    zeros[i] ~ dpois(zeros.mean[i])
    
    ##zero trick
    zeros.mean[i] <- -l[i] + C
    
    ##customize log-likelihood of the weibull distribution
    l[i] <- log(ifelse(is.observed[i],
            h[i] * S[i], 
            S[i]))
            
    ##constructing hazard and Survival function
    beta.form[i] <- beta[1] * age[i] + 
                    beta[2] * sex[i] +
                    beta[3] * body[i] +
                    beta[4] * tcic[i] + 
                    beta[5] * smoke[i]
    
    delta.form[i] <- delta[1] * age[i] * trt[i] + 
                     delta[2] * sex[i] * trt[i] +
                     delta[3] * body[i] * trt[i] +
                     delta[4] * tcic[i] * trt[i] + 
                     delta[5] * smoke[i] * trt[i] + 
                     delta[6] * smoke[i] * tcic[i] * trt[i] + 
                     delta[7] * age[i] * body[i] * trt[i]
    
    ##when treatment = T and Time_i >= tau       
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau + tau * exp(alpha + delta.form[i])))
                              - mu - alpha - delta.form[i] - beta.form[i])
    S.tau[i] <- pow(1 + exp(epsilon.form.tau[i]), -1)
    h.tau[i] <- b / (abs(time[i] - tau + tau * exp(alpha + delta.form[i])))
                    * pow(1 + exp(-epsilon.form.tau[i]), -1)
    
    ##otherwise
    epsilon.form.ll[i] <- b * (log(time[i]) - mu - beta.form[i])
    S.ll[i] <- pow(1 + exp(epsilon.form.ll[i]), -1)
    h.ll[i] <- b/time[i] * pow(1 + exp(-epsilon.form.ll[i]), -1)
    
    S[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, S.tau[i], S.ll[i]),
                   S.ll[i])
    h[i] <- ifelse(trt[i],
                   ifelse(time[i] > tau, h.tau[i], h.ll[i]),
                   h.ll[i])
  }
  
  ##priors for betas
  alpha ~ dnorm(0.35553, 75)
  beta[1] ~ dnorm(-0.05423850, 250)
  beta[2] ~ dnorm(0.16185739, 100)
  beta[3] ~ dnorm(-0.60899904, 100)
  beta[4] ~ dnorm(0.09379640, 100)
  beta[5] ~ dnorm(-0.06814758, 50)
  
  for(i in 1:7){
    delta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the log-logistic distribution
  mu ~ dnorm(2.60503197, 10)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  theta <- -mu/sigma
  tau ~ dnorm(2.01429825, 1)I(0,5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.cusLL.delta2 = file.path(dir, 'cusLL.delta2.txt')
writeLines(modelstring.cusLL.delta2, con=modelfile.cusLL.delta2)

set.seed(123)
auto.ll.delta3 = autorun.jags(model = modelfile.cusLL.delta2, data = datlist.tau, 
                              monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                              n.chains = 3, startburnin = 5000, startsample = 6000, 
                              thin = 1, method = 'parallel')
auto.ll.delta2
auto.ll.delta3




########################################################################
########################################################################
########################################################################
##plots used for 6-page black and white shorter paper

ggsurv(sf, surv.col = c('black', 'gray')) + geom_vline(xintercept = 2.7) + 
  coord_cartesian(xlim = c(0, 30)) + 
  theme(plot.title = element_text(hjust = 0.6, face = 'bold', size = 20)) + 
  theme_bw()


plotbw.xy1 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=age)) + 
  scale_color_gradient(low = "black", high = "gray90") + 
  labs(x = 'PRT', y = 'TE', color = 'Age') +
  theme_bw()
##by sex
plotbw.xy2 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=sex)) + 
  scale_color_grey() + 
  labs(x = 'PRT', y = 'TE', color = 'Sex') +
  theme_bw()
##by body
plotbw.xy3 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=body)) + 
  scale_color_grey() + 
  labs(x = 'PRT', y = 'TE', color = 'Body') +
  theme_bw()
##by tcic
plotbw.xy4 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=tcic)) + 
  scale_color_grey() + 
  labs(x = 'PRT', y = 'TE', color = 'TCIC') +
  theme_bw()
##by smoke
plotbw.xy5 = ggplot(Xplot[nonzero,], aes(x=prob.morethantau, y=delta.treat)) + 
  geom_point(aes(color=smoke)) + 
  scale_color_grey() + 
  labs(x = 'PRT', y = 'TE', color = 'Smoke') +
  theme_bw()
##cluster analysis
plotbw.ca = ggplot(as.data.frame(cluster.plot), aes(x = `PRT`, 
                                                    y = `TE`, 
                                                    color = Cluster5)) +
  geom_point(size = 1.5) +
  scale_color_grey(start = 0, end = 0.9) +  # adjust the range for better contrast
  labs(x = "PRT", y = "TE", color = '5 Clusters') +
  theme_bw()
plotbw.xy = ggarrange(plotbw.xy1, plotbw.xy2, plotbw.xy3, plotbw.xy4, 
                      plotbw.xy5, plotbw.ca, ncol=3, nrow=2)
annotate_figure(plotbw.xy)


pheatmap(t(cluster.std),
         color = gray.colors(100, start = 0, end = 0.95),
         cluster_cols = hc,
         cluster_rows = FALSE,
         clustering_method = 'average',
         cutree_cols = 5,
         treeheight_col = 100,
         annotation_names_row = 1
)





