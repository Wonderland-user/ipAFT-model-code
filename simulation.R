library(rjags)
library(DescTools)

##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
alpha0 = 1.5; tau0 = 1; beta10 = 2; beta20 = -1.8


##treatment Z ~ Bin(n, 0.5)
##x1 ~ N(0.6, 0.4)
##x2 ~ LN(-0.8, 0.6)
##epsilon ~ Gumbel(0, 1)
##set n = 100
n = 400
Z = rbernoulli(n=n, p=0.5)
x1 = rnorm(n = n, mean = 0.6, sd = sqrt(0.4))
x2 = rlnorm(n = n, meanlog = -0.8, sdlog = sqrt(0.6))
epsilon = log(rexp(n))


##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
tauT = T0

##now determine which patient has treatment effect and change their T
for(i in 1:n){
  if(Z[i]){
    if(T0[i] > tau0){
      ##patients with treatment and T > tau0
      tauT[i] = tau0 + exp(alpha0) * (T0[i] - tau0)
    }
  }
}

##apply censoring. 25% of the patients will be censored
for(i in 1:(n/4)){
  tauT[i] = runif(n = 1, min = 0, max = tauT[i])
}
is.observed.tau = c(rep(FALSE, times = n/4), rep(TRUE, times = 3*n/4))

##set a trial length maximum at 500
for(i in 1:n){
  if(tauT[i] > 700){
    tauT[i] = 700
    is.observed.tau[i] = FALSE
  }
}


##now check if tau model work
modelstring.tauModel = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    ##when treatment = T and Time_i > tau        
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
  
  ##priors
  alpha ~ dnorm(0, 0.001)
  for(i in 1:2){
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
modelfile.tauModel = file.path(dir, 'tauModel.txt')
writeLines(modelstring.tauModel, con=modelfile.tauModel)
datlist.tau = list(N = n,
               time = tauT,
               is.observed = is.observed.tau,
               trt = Z,
               x1 = x1,
               x2 = x2,
               pi = pi)

tauModel = jags.parallel(data = datlist.tau,
                         parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                         n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                         n.thin = 1, model.file = modelfile.tauModel)
tauModel$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(tauModel)


########################################################################################
########################################################################################
##Eta model
length(etaT[etaT!=T0])
##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
alpha0 = 1; tau0 = 3; beta10 = -0.5; beta20 = 0.7; mu0=2
##set eta1 = 1.5 and eta2 = 0.2
eta10 = 0.3; eta20 = -0.3

##treatment Z ~ Bin(n, 0.5)
##x1 ~ N(0.6, 0.4)
##x2 ~ LN(-0.8, 0.6)
##epsilon ~ Gumbel(0, 1)
##set n = 100
n = 800
Z = rbernoulli(n=n, p=0.5)
x1 = rnorm(n = n, mean = 0.6, sd = sqrt(0.4))
x2 = rlnorm(n = n, meanlog = -0.8, sdlog = sqrt(0.6))
x2 = (x2-mean(x2))/sd(x2) ## standardize x2
#x2 = rbernoulli(n=n, p=0.5)

epsilon = log(rexp(n))

##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
T0 = exp(mu0 + beta10 * x1 + beta20 * x2 + epsilon)

##now determine which patient has treatment effect and change their T
tau0.eta = tau0^exp(eta10 * x1 + eta20 * x2)
etaT = T0
for(i in 1:n){
  if(Z[i]){
    if(T0[i] > tau0.eta[i]){
      ##patients with treatment and T > tau0
      etaT[i] = tau0.eta[i] + exp(alpha0) * (T0[i] - tau0.eta[i])
    }
  }
}

##apply censoring. 25% of the patients will be censored
for(i in 1:(n/4)){
  etaT[i] = runif(n = 1, min = 0, max = etaT[i])
}
is.observed.eta = c(rep(FALSE, times = n/4), rep(TRUE, times = 3*n/4))

##set a trial length maximum at 700
for(i in 1:n){
  if(etaT[i] > 700){
    etaT[i] = 700
    is.observed.eta[i] = FALSE
  }
}
Treatment = Z

sf = survfit(Surv(etaT, is.observed.eta) ~ Treatment)
ggsurv(sf) + coord_cartesian(xlim = c(0, 50)) + 
  labs(title = 'Kaplan-Meier Curve of a Simulated Dataset') + 
  theme(plot.title = element_text(hjust = 0.6, face = 'bold', size = 20))


##test eta model
modelstring.etaModel = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    eta.form[i] <- pow(tau, exp(eta[1] * x1[i] + eta[2] * x2[i]))
    
    ##when treatment = T and Time_i >= exp(tau + eta*Z)       
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
  
  alpha ~ dnorm(1, 5)
  beta[1] ~ dnorm(-0.5, 5)
  beta[2] ~ dnorm(0.7, 5)
  #alpha ~ dnorm(0, 0.01)
  #beta[1] ~ dnorm(0, 0.01)
  #beta[2] ~ dnorm(0, 0.01)
  
  for(i in 1:2){
    eta[i] ~ dnorm(0, 1/2)
  }

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.01)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(3, 1)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
#dir = 'C:/Users/xingz/OneDrive/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.etaModel = file.path(dir, 'etaModel.txt')
writeLines(modelstring.etaModel, con=modelfile.etaModel)

datlist.eta = list(N = n,
                   time = etaT,
                   is.observed = is.observed.eta,
                   trt = Z,
                   x1 = x1,
                   x2 = x2,
                   pi = pi)

etaModel = jags.parallel(data = datlist.eta, 
                         parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'eta', 'alpha'),
                         n.chains = 5, n.burnin = 7000, n.iter = 9000, 
                         n.thin = 1, model.file = modelfile.etaModel)
etaModel$BUGSoutput[10]$summary[,c(1,3,7)]

datlist.eta2 = list(N = n,
                    time = etaT,
                    is.observed = ifelse(is.observed.eta, 1, 0),
                    trt = ifelse(Z, 1, 0),
                    x1 = x1,
                    x2 = as.integer(x2),
                    pi = pi)
etaM800 = autorun.jags(model = modelfile.etaModel, data = datlist.eta2, 
                    monitor = c('alpha', 'beta', 'mu', 'sigma', 'tau', 'eta'), 
                    n.chains = 3, startburnin = 5000, startsample = 6000, 
                    thin = 1, method = 'parallel')
etaM
etaM800

etaM800.2 = autorun.jags(model = modelfile.etaModel, data = datlist.eta2, 
                         monitor = c('alpha', 'beta', 'mu', 'sigma', 'tau', 'eta'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 10, method = 'parallel')
etaM800.2
#traceplot(etaModel)

#This model cannot estimate eta1 and eta2.




########################################################################################
########################################################################################
##Delta Model

##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
##set delta1 = 4 and delta2 = 0.2
delta10 = -1.5; delta20 = 1.5

##treatment Z ~ Bin(n, 0.5)
##x1 ~ N(0.6, 0.4)
##x2 ~ LN(-0.8, 0.6)
##epsilon ~ Gumbel(0, 1)
##set n = 100
n = 800
Z = rbernoulli(n=n, p=0.5)
x1 = rnorm(n = n, mean = 0.6, sd = sqrt(0.4))
x2 = rlnorm(n = n, meanlog = -0.8, sdlog = sqrt(0.6))
epsilon = log(rexp(n))

##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)

##now determine which patient has treatment effect and change their T
alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
deltaT = T0
for(i in 1:n){
  if(Z[i]){
    if(T0[i] > tau0){
      ##patients with treatment and T > tau0
      deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
    }
  }
}

##apply censoring. 25% of the patients will be censored
for(i in 1:(n/4)){
  deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
}
is.observed.delta = c(rep(FALSE, times = n/4), rep(TRUE, times = 3*n/4))

##set a trial length maximum at 500
for(i in 1:n){
  if(deltaT[i] > 700){
    deltaT[i] = 700
    is.observed.delta[i] = FALSE
  }
}


modelstring.deltaModel = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  #alpha ~ dnorm(1.506882, 5)
  #beta[1] ~ dnorm(2.004770, 5)
  #beta[2] ~ dnorm(-1.807582, 5)
  #alpha ~ dnorm(0.4469, 1)
  #beta[1] ~ dnorm(1.9624, 5)
  #beta[2] ~ dnorm(-1.6884, 5)
  alpha ~ dnorm(0, 0.001)
  beta[1] ~ dnorm(0, 0.001)
  beta[2] ~ dnorm(0, 0.001)
   
  for(i in 1:2){
    delta[i] ~ dnorm(0, 0.001)
  }

  ##prior for mu and sigma of the weibull distribution
  #mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  #tau ~ dnorm(2.553775, 10)I(0,)
  
  #mu ~ dnorm(0.16188, 0.001)
  mu ~ dnorm(0, 0.001)
  tau ~ dunif(0.5, 5)
  #tau ~ dnorm(2.2174, 1)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.deltaModel = file.path(dir, 'deltaModel.txt')
writeLines(modelstring.deltaModel, con=modelfile.deltaModel)

datlist.delta = list(N = n,
                   time = deltaT,
                   is.observed = is.observed.delta,
                   trt = Z,
                   x1 = x1,
                   x2 = x2,
                   pi = pi)

deltaModel = jags.parallel(data = datlist.delta, 
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'delta', 'alpha'),
                           n.chains = 3, n.burnin = 6000, n.iter = 8000, 
                           n.thin = 1, model.file = modelfile.deltaModel)
deltaModel$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(deltaModel)

##use tauModel to restrict priors
##prior = cbind(tauModel$BUGSoutput[10]$summary[,c(1,2)], 1/(tauModel$BUGSoutput[10]$summary[,2])^2)
##prior


datlist.delta = list(N = n,
                     time = deltaT,
                     is.observed = as.integer(is.observed.delta),
                     trt = as.integer(Z),
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
ipAFT_vag = autorun.jags(model = modelfile.deltaModel, data = datlist.delta, 
                         monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 1, method = 'parallel')

pAFT = autorun.jags(model = modelfile.tauModel, data = datlist.delta, 
                    monitor = c('alpha', 'beta', 'mu', 'sigma', 'tau'), 
                    n.chains = 3, startburnin = 5000, startsample = 6000, 
                    thin = 1, method = 'parallel')
ipAFT_inf = autorun.jags(model = modelfile.deltaModel, data = datlist.delta, 
                         monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 1, method = 'parallel')
#res = ipAFT_inf
time_vag = system.time(autorun.jags(model = modelfile.deltaModel, data = datlist.delta, 
                         monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 1, method = 'parallel'))
time_inf = system.time(autorun.jags(model = modelfile.tauModel, data = datlist.delta, 
                         monitor = c('alpha', 'beta', 'delta', 'mu', 'sigma', 'tau'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 1, method = 'parallel'))

########################################################################################
########################################################################################
##Convergence Check
library(coda)
delta.mcmc = as.mcmc(deltaModel)
gelman.diag(delta.mcmc)


eta.mcmc = as.mcmc(etaModel)
gelman.diag(eta.mcmc)
#gelman.plot(eta.mcmc)
########################################################################################
########################################################################################
##KM plot
library(survival)
library(ggplot2)
library(GGally)

set.seed(123)
##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 2; delta20 = -2
n = 400; W = rbernoulli(n=n, p=0.5); epsilon = log(rexp(n))
x1 = rnorm(n = n, mean = 0.6, sd = sqrt(0.4))
x2 = rlnorm(n = n, meanlog = -0.8, sdlog = sqrt(0.6))

T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
alpha0.plot = alpha0 + delta10 * x1 + delta20 * x2
plotT = T0
for(i in 1:n){
  if(W[i]){
    if(T0[i] > tau0){
      ##patients with treatment and T > tau0
      plotT[i] = tau0 + exp(alpha0.plot[i]) * (T0[i] - tau0)
    }
  }
}


for(i in 1:(n/4)){
  plotT[i] = runif(n = 1, min = 0, max = plotT[i])
}
is.observed.plot = c(rep(FALSE, times = n/4), rep(TRUE, times = 3*n/4))

for(i in 1:n){
  if(plotT[i] > 700){
    plotT[i] = 700
    is.observed.plot[i] = FALSE
  }
}

sf = survfit(Surv(plotT, is.observed.plot) ~ W)
ggsurv(sf) + coord_cartesian(xlim = c(0, 10))

datlist.plot = list(N = n,
                    time = plotT,
                    is.observed = is.observed.plot,
                    trt = Z,
                    x1 = x1,
                    x2 = x2,
                    pi = pi)
plotModel = jags.parallel(data = datlist.plot,
                          parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                          n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                          n.thin = 1, model.file = modelfile.tauModel)
plotModel$BUGSoutput[10]$summary[,c(1,3,7)]

##
newT = plotT[plotT >= 2.51] - 2.51
newCen = is.observed.plot[plotT >= 2.51]
newZ = Z[plotT >= 2.51]

newsf = survfit(Surv(newT, newCen) ~ newZ)
ggsurv(newsf) + coord_cartesian(xlim = c(0, 100))

newcox = coxph(Surv(newT, newCen) ~ newZ)
newcox

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

##Simulation study
##n=400; 100 times each

##Tau Model
set.seed(111)
npatient = 400
nsim = 100
tau.sim = c()
tau.mu = rep(0, times = nsim)
trueVal.tau = c(alpha0, beta10, beta20, tau0)
tau.cover = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.cover) = c('alpha', 'beta1', 'beta2', 'tau')
tau.order = c(1, 2, 3, 7)
tau.significance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.significance) = c('alpha', 'beta1', 'beta2', 'tau')
tau.wrongsignificance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('tau current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  tauT = T0
  
  ##now determine which patient has treatment effect and change their T
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        tauT[i] = tau0 + exp(alpha0) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    tauT[i] = runif(n = 1, min = 0, max = tauT[i])
  }
  is.observed.tau = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(tauT[i] > 700){
      tauT[i] = 700
      is.observed.tau[i] = FALSE
    }
  }
  
  datlist.tau = list(N = npatient,
                     time = tauT,
                     is.observed = is.observed.tau,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  tauModel = jags.parallel(data = datlist.tau,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                           n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                           n.thin = 1, model.file = modelfile.tauModel)
  tau.sim = cbind(tau.sim, tauModel$BUGSoutput[10]$summary[c(1:3,7),1])
  
  ##check if this run covers the true coef
  for(i in 1:4){
    if(trueVal.tau[i] >= tauModel$BUGSoutput[10]$summary[tau.order[i], 3]){
      if(trueVal.tau[i] <= tauModel$BUGSoutput[10]$summary[tau.order[i], 7]){
        tau.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(tauModel$BUGSoutput[10]$summary[1, 3] >= 0){tau.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(tauModel$BUGSoutput[10]$summary[1, 7] <= 0){tau.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(tauModel$BUGSoutput[10]$summary[2, 3] >= 0){tau.significance[2, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[2, 7] <= 0){tau.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(tauModel$BUGSoutput[10]$summary[3, 7] <= 0){tau.significance[3, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[3, 3] >= 0){tau.wrongsignificance[3, isim] = TRUE}
  
  ##tau
  if(tauModel$BUGSoutput[10]$summary[7, 3] >= 0){tau.significance[4, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[7, 7] <= 0){tau.wrongsignificance[4, isim] = TRUE}
  
  ##record mu for setting the prior of delta models
  tau.mu[isim] = tauModel$BUGSoutput[10]$summary[5,1]
}
##mean of the simulation of tau model
(tau.mean = apply(tau.sim, FUN = mean, 1))

##mean of the bias of the simulated tau model
trueVal.tau = c(alpha0, beta10, beta20, tau0)
(tau.bias = apply(tau.sim - trueVal.tau, FUN = mean, 1))

mean(tau.mu)

##result
tau400 = cbind(trueVal.tau, tau.mean, tau.bias)
tau.ci = rbind(as.numeric(t.test(tau.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[4,], conf.level = 0.95)$conf.int))
tau400 = cbind(tau400, tau.ci, apply(tau.cover, FUN = mean, 1),
               apply(tau.significance, FUN = mean, 1),
               apply(tau.wrongsignificance, FUN = mean, 1))
tau.sim400 = tau.sim
colnames(tau400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
tau400


########################################################################################
########################################################################################
##simulation for delta model
set.seed(111)
npatient = 400
nsim = 100
delta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                     time = deltaT,
                     is.observed = is.observed.delta,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  deltaModel = jags.parallel(data = datlist.delta,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                           n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                           n.thin = 1, model.file = modelfile.deltaModel)
  delta.sim = cbind(delta.sim, deltaModel$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= deltaModel$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= deltaModel$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 3] >= 0){delta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 7] <= 0){delta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(deltaModel$BUGSoutput[10]$summary[2, 3] >= 0){delta.significance[2, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[2, 7] <= 0){delta.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(deltaModel$BUGSoutput[10]$summary[3, 7] <= 0){delta.significance[3, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[3, 3] >= 0){delta.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(deltaModel$BUGSoutput[10]$summary[4, 3] >= 0){delta.significance[4, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[4, 7] <= 0){delta.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(deltaModel$BUGSoutput[10]$summary[5, 7] <= 0){delta.significance[5, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[5, 3] >= 0){delta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(deltaModel$BUGSoutput[10]$summary[9, 3] >= 0){delta.significance[6, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[9, 7] <= 0){delta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.mean = apply(delta.sim, FUN = mean, 1))
(delta.bias = apply(delta.sim - trueVal.delta, FUN = mean, 1))

delta400 = cbind(trueVal.delta, delta.mean, delta.bias)
delta.ci = rbind(as.numeric(t.test(delta.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(delta.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(delta.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(delta.sim[4,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(delta.sim[5,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(delta.sim[6,], conf.level = 0.95)$conf.int))
delta400 = cbind(delta400, delta.ci, apply(delta.cover, FUN = mean, 1),
                 apply(delta.significance, FUN = mean, 1),
                 apply(delta.wrongsignificance, FUN = mean, 1))
colnames(delta400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta400
##hard to converge
#deltaModel$BUGSoutput[10]$summary[,c(1,3,7)]
#traceplot(deltaModel)


##
##set tau beta1 beta2 prior to tau model. Also look for potential prior for delta1 and delta2
##
##divide exp(alpha + (beta + delta)X^t)*(T_0i - tau) by exp(beta*X^t)*(T_0j - tau).
##the result can be used to estimate delta. (alpha, beta, T_0i, T_0j, tau are estimated in tau model)

##get 95CI for tau400 and delta400.

##number of runs true coef is covered in the 95% confint




########################################################################################
########################################################################################
##simulation for eta model
set.seed(111)
npatient = 400
nsim = 100
eta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; eta10 = -1; eta20 = 1
trueVal.eta = c(alpha0, beta10, beta20, eta10, eta20, tau0)
eta.cover = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.cover) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.order = c(1:3, 5:6, 9)
eta.significance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.significance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('eta current run', isim))
  
  ##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  ##set eta1 = 1.5 and eta2 = 0.2
  eta10 = -1; eta20 = 1
  
  ##treatment Z ~ Bin(n, 0.5)
  ##x1 ~ N(0.6, 0.4)
  ##x2 ~ LN(-0.8, 0.6)
  ##epsilon ~ Gumbel(0, 1)
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  
  ##now determine which patient has treatment effect and change their T
  tau0.eta = tau0 + eta10 * x1 + eta20 * x2
  etaT = T0
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0.eta[i]){
        ##patients with treatment and T > tau0
        etaT[i] = tau0.eta[i] + exp(alpha0) * (T0[i] - tau0.eta[i])
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    etaT[i] = runif(n = 1, min = 0, max = etaT[i])
  }
  is.observed.eta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(etaT[i] > 700){
      etaT[i] = 700
      is.observed.eta[i] = FALSE
    }
  }
  
  datlist.eta = list(N = npatient,
                       time = etaT,
                       is.observed = is.observed.eta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  etaModel = jags.parallel(data = datlist.eta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'),
                             n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                             n.thin = 1, model.file = modelfile.etaModel)
  eta.sim = cbind(eta.sim, etaModel$BUGSoutput[10]$summary[c(1:3, 5:6,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.eta)){
    if(trueVal.eta[i] >= etaModel$BUGSoutput[10]$summary[eta.order[i], 3]){
      if(trueVal.eta[i] <= etaModel$BUGSoutput[10]$summary[eta.order[i], 7]){
        eta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(etaModel$BUGSoutput[10]$summary[1, 3] >= 0){eta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(etaModel$BUGSoutput[10]$summary[1, 7] <= 0){eta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(etaModel$BUGSoutput[10]$summary[2, 3] >= 0){eta.significance[2, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[2, 7] <= 0){eta.wrongsignificance[2, isim] = TRUE}

  ##beta2
  if(etaModel$BUGSoutput[10]$summary[3, 7] <= 0){eta.significance[3, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[3, 3] >= 0){eta.wrongsignificance[3, isim] = TRUE}
  
  ##eta1
  if(etaModel$BUGSoutput[10]$summary[5, 7] <= 0){eta.significance[4, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[5, 3] >= 0){eta.wrongsignificance[4, isim] = TRUE}
  
  ##eta2
  if(etaModel$BUGSoutput[10]$summary[6, 3] >= 0){eta.significance[5, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[6, 7] <= 0){eta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(etaModel$BUGSoutput[10]$summary[9, 3] >= 0){eta.significance[6, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[9, 7] <= 0){eta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated eta model
trueVal.eta = c(alpha0, beta10, beta20, eta10, eta20, tau0)
(eta.mean = apply(eta.sim, FUN = mean, 1))
(eta.bias = apply(eta.sim - trueVal.eta, FUN = mean, 1))

eta400 = cbind(trueVal.eta, eta.mean, eta.bias)
eta.ci = rbind(as.numeric(t.test(eta.sim[1,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(eta.sim[2,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(eta.sim[3,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(eta.sim[4,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(eta.sim[5,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(eta.sim[6,], conf.level = 0.95)$conf.int))
eta400 = cbind(eta400, eta.ci, apply(eta.cover, FUN = mean, 1), 
               apply(eta.significance, FUN = mean, 1),
               apply(eta.wrongsignificance, FUN = mean, 1))
colnames(eta400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
eta400















########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

##Simulation study
##800

##Tau Model
set.seed(123)
npatient = 800
nsim = 1000
tau.sim = c()
tau.mu = rep(0, times = nsim)
trueVal.tau = c(alpha0, beta10, beta20, tau0)
tau.cover = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.cover) = c('alpha', 'beta1', 'beta2', 'tau')
tau.order = c(1, 2, 3, 7)
tau.significance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.significance) = c('alpha', 'beta1', 'beta2', 'tau')
tau.wrongsignificance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('tau current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  tauT = T0
  
  ##now determine which patient has treatment effect and change their T
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        tauT[i] = tau0 + exp(alpha0) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    tauT[i] = runif(n = 1, min = 0, max = tauT[i])
  }
  is.observed.tau = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(tauT[i] > 700){
      tauT[i] = 700
      is.observed.tau[i] = FALSE
    }
  }
  
  datlist.tau = list(N = npatient,
                     time = tauT,
                     is.observed = is.observed.tau,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  tauModel = jags.parallel(data = datlist.tau,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                           n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                           n.thin = 1, model.file = modelfile.tauModel)
  tau.sim = cbind(tau.sim, tauModel$BUGSoutput[10]$summary[c(1:3,7),1])
  
  ##check if this run covers the true coef
  for(i in 1:4){
    if(trueVal.tau[i] >= tauModel$BUGSoutput[10]$summary[tau.order[i], 3]){
      if(trueVal.tau[i] <= tauModel$BUGSoutput[10]$summary[tau.order[i], 7]){
        tau.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(tauModel$BUGSoutput[10]$summary[1, 3] >= 0){tau.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(tauModel$BUGSoutput[10]$summary[1, 7] <= 0){tau.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(tauModel$BUGSoutput[10]$summary[2, 3] >= 0){tau.significance[2, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[2, 7] <= 0){tau.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(tauModel$BUGSoutput[10]$summary[3, 7] <= 0){tau.significance[3, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[3, 3] >= 0){tau.wrongsignificance[3, isim] = TRUE}
  
  ##tau
  if(tauModel$BUGSoutput[10]$summary[7, 3] >= 0){tau.significance[4, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[7, 7] <= 0){tau.wrongsignificance[4, isim] = TRUE}
  
  ##record mu for setting the prior of delta models
  tau.mu[isim] = tauModel$BUGSoutput[10]$summary[5,1]
}
##mean of the simulation of tau model
(tau.mean = apply(tau.sim, FUN = mean, 1))

##mean of the bias of the simulated tau model
trueVal.tau = c(alpha0, beta10, beta20, tau0)
(tau.bias = apply(tau.sim - trueVal.tau, FUN = mean, 1))

mean(tau.mu)

##result
tau800 = cbind(trueVal.tau, tau.mean, tau.bias)
tau.ci = rbind(as.numeric(t.test(tau.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[4,], conf.level = 0.95)$conf.int))
tau800 = cbind(tau800, tau.ci, apply(tau.cover, FUN = mean, 1),
               apply(tau.significance, FUN = mean, 1),
               apply(tau.wrongsignificance, FUN = mean, 1))
tau.sim800 = tau.sim
colnames(tau800) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
tau800


########################################################################################
########################################################################################
##simulation for delta model
set.seed(111)
npatient = 800
nsim = 1000
delta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  deltaModel = jags.parallel(data = datlist.delta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                             n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                             n.thin = 1, model.file = modelfile.deltaModel)
  delta.sim = cbind(delta.sim, deltaModel$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= deltaModel$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= deltaModel$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 3] >= 0){delta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 7] <= 0){delta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(deltaModel$BUGSoutput[10]$summary[2, 3] >= 0){delta.significance[2, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[2, 7] <= 0){delta.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(deltaModel$BUGSoutput[10]$summary[3, 7] <= 0){delta.significance[3, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[3, 3] >= 0){delta.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(deltaModel$BUGSoutput[10]$summary[4, 3] >= 0){delta.significance[4, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[4, 7] <= 0){delta.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(deltaModel$BUGSoutput[10]$summary[5, 7] <= 0){delta.significance[5, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[5, 3] >= 0){delta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(deltaModel$BUGSoutput[10]$summary[9, 3] >= 0){delta.significance[6, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[9, 7] <= 0){delta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.mean = apply(delta.sim, FUN = mean, 1))
(delta.bias = apply(delta.sim - trueVal.delta, FUN = mean, 1))

delta800 = cbind(trueVal.delta, delta.mean, delta.bias)
delta.ci = rbind(as.numeric(t.test(delta.sim[1,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[2,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[3,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[4,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[5,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[6,], conf.level = 0.95)$conf.int))
delta800 = cbind(delta800, delta.ci, apply(delta.cover, FUN = mean, 1),
                 apply(delta.significance, FUN = mean, 1),
                 apply(delta.wrongsignificance, FUN = mean, 1))
colnames(delta800) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta800





########################################################################################
########################################################################################
##simulation for eta model
set.seed(111)
npatient = 800
nsim = 100
eta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; eta10 = -1; eta20 = 1
trueVal.eta = c(alpha0, beta10, beta20, eta10, eta20, tau0)
eta.cover = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.cover) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.order = c(1:3, 5:6, 9)
eta.significance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.significance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('eta current run', isim))
  
  ##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  ##set eta1 = 1.5 and eta2 = 0.2
  eta10 = -1; eta20 = 1
  
  ##treatment Z ~ Bin(n, 0.5)
  ##x1 ~ N(0.6, 0.4)
  ##x2 ~ LN(-0.8, 0.6)
  ##epsilon ~ Gumbel(0, 1)
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  
  ##now determine which patient has treatment effect and change their T
  tau0.eta = tau0 + eta10 * x1 + eta20 * x2
  etaT = T0
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0.eta[i]){
        ##patients with treatment and T > tau0
        etaT[i] = tau0.eta[i] + exp(alpha0) * (T0[i] - tau0.eta[i])
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    etaT[i] = runif(n = 1, min = 0, max = etaT[i])
  }
  is.observed.eta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(etaT[i] > 700){
      etaT[i] = 700
      is.observed.eta[i] = FALSE
    }
  }
  
  datlist.eta = list(N = npatient,
                     time = etaT,
                     is.observed = is.observed.eta,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  etaModel = jags.parallel(data = datlist.eta,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'),
                           n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                           n.thin = 1, model.file = modelfile.etaModel)
  eta.sim = cbind(eta.sim, etaModel$BUGSoutput[10]$summary[c(1:3, 5:6,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.eta)){
    if(trueVal.eta[i] >= etaModel$BUGSoutput[10]$summary[eta.order[i], 3]){
      if(trueVal.eta[i] <= etaModel$BUGSoutput[10]$summary[eta.order[i], 7]){
        eta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(etaModel$BUGSoutput[10]$summary[1, 3] >= 0){eta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(etaModel$BUGSoutput[10]$summary[1, 7] <= 0){eta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(etaModel$BUGSoutput[10]$summary[2, 3] >= 0){eta.significance[2, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[2, 7] <= 0){eta.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(etaModel$BUGSoutput[10]$summary[3, 7] <= 0){eta.significance[3, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[3, 3] >= 0){eta.wrongsignificance[3, isim] = TRUE}
  
  ##eta1
  if(etaModel$BUGSoutput[10]$summary[5, 7] <= 0){eta.significance[4, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[5, 3] >= 0){eta.wrongsignificance[4, isim] = TRUE}
  
  ##eta2
  if(etaModel$BUGSoutput[10]$summary[6, 3] >= 0){eta.significance[5, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[6, 7] <= 0){eta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(etaModel$BUGSoutput[10]$summary[9, 3] >= 0){eta.significance[6, isim] = TRUE}
  if(etaModel$BUGSoutput[10]$summary[9, 7] <= 0){eta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated eta model
trueVal.eta = c(alpha0, beta10, beta20, eta10, eta20, tau0)
(eta.mean = apply(eta.sim, FUN = mean, 1))
(eta.bias = apply(eta.sim - trueVal.eta, FUN = mean, 1))

eta800 = cbind(trueVal.eta, eta.mean, eta.bias)
eta.ci = rbind(as.numeric(t.test(eta.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.sim[4,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.sim[5,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.sim[6,], conf.level = 0.95)$conf.int))
eta800 = cbind(eta800, eta.ci, apply(eta.cover, FUN = mean, 1), 
               apply(eta.significance, FUN = mean, 1),
               apply(eta.wrongsignificance, FUN = mean, 1))
colnames(eta800) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
eta800























########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

##Simulation study
##n=1600; 100 times each

##Tau Model
set.seed(111)
npatient = 1600
nsim = 100
tau.sim = c()
tau.mu = rep(0, times = nsim)
trueVal.tau = c(alpha0, beta10, beta20, tau0)
tau.cover = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.cover) = c('alpha', 'beta1', 'beta2', 'tau')
tau.order = c(1, 2, 3, 7)
tau.significance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.significance) = c('alpha', 'beta1', 'beta2', 'tau')
tau.wrongsignificance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('tau current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  tauT = T0
  
  ##now determine which patient has treatment effect and change their T
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        tauT[i] = tau0 + exp(alpha0) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    tauT[i] = runif(n = 1, min = 0, max = tauT[i])
  }
  is.observed.tau = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(tauT[i] > 700){
      tauT[i] = 700
      is.observed.tau[i] = FALSE
    }
  }
  
  datlist.tau = list(N = npatient,
                     time = tauT,
                     is.observed = is.observed.tau,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  tauModel = jags.parallel(data = datlist.tau,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                           n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                           n.thin = 1, model.file = modelfile.tauModel)
  tau.sim = cbind(tau.sim, tauModel$BUGSoutput[10]$summary[c(1:3,7),1])
  
  ##check if this run covers the true coef
  for(i in 1:4){
    if(trueVal.tau[i] >= tauModel$BUGSoutput[10]$summary[tau.order[i], 3]){
      if(trueVal.tau[i] <= tauModel$BUGSoutput[10]$summary[tau.order[i], 7]){
        tau.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(tauModel$BUGSoutput[10]$summary[1, 3] >= 0){tau.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(tauModel$BUGSoutput[10]$summary[1, 7] <= 0){tau.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(tauModel$BUGSoutput[10]$summary[2, 3] >= 0){tau.significance[2, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[2, 7] <= 0){tau.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(tauModel$BUGSoutput[10]$summary[3, 7] <= 0){tau.significance[3, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[3, 3] >= 0){tau.wrongsignificance[3, isim] = TRUE}
  
  ##tau
  if(tauModel$BUGSoutput[10]$summary[7, 3] >= 0){tau.significance[4, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[7, 7] <= 0){tau.wrongsignificance[4, isim] = TRUE}
  
  ##record mu for setting the prior of delta models
  tau.mu[isim] = tauModel$BUGSoutput[10]$summary[5,1]
}
##mean of the simulation of tau model
(tau.mean = apply(tau.sim, FUN = mean, 1))

##mean of the bias of the simulated tau model
trueVal.tau = c(alpha0, beta10, beta20, tau0)
(tau.bias = apply(tau.sim - trueVal.tau, FUN = mean, 1))

mean(tau.mu)

##result
tau1600 = cbind(trueVal.tau, tau.mean, tau.bias)
tau.ci = rbind(as.numeric(t.test(tau.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[4,], conf.level = 0.95)$conf.int))
tau1600 = cbind(tau1600, tau.ci, apply(tau.cover, FUN = mean, 1),
               apply(tau.significance, FUN = mean, 1),
               apply(tau.wrongsignificance, FUN = mean, 1))
tau.sim1600 = tau.sim
colnames(tau1600) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
tau1600













########################################################################################
########################################################################################
##simulation for eeta model
set.seed(111)
npatient = 400
nsim = 5
eeta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; eeta10 = -2; eeta20 = 2
trueVal.eeta = c(alpha0, beta10, beta20, eeta10, eeta20, tau0)
eeta.cover = matrix(FALSE, nrow = length(trueVal.eeta), ncol = nsim)
rownames(eeta.cover) = c('alpha', 'beta1', 'beta2', 'eeta1', 'eeta2', 'tau')
eeta.order = c(1:3, 5:6, 9)
eeta.significance = matrix(FALSE, nrow = length(trueVal.eeta), ncol = nsim)
rownames(eeta.significance) = c('alpha', 'beta1', 'beta2', 'eeta1', 'eeta2', 'tau')
eeta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.eeta), ncol = nsim)
rownames(eeta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'eeta1', 'eeta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('eeta current run', isim))
  
  ##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  eeta10 = -2; eeta20 = 2
  
  ##treatment Z ~ Bin(n, 0.5)
  ##x1 ~ N(0.6, 0.4)
  ##x2 ~ LN(-0.8, 0.6)
  ##epsilon ~ Gumbel(0, 1)
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  
  ##now determine which patient has treatment effect and change their T
  tau0.eeta = tau0 + eeta10 * x1 + eeta20 * x2
  eetaT = T0
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0.eeta[i]){
        ##patients with treatment and T > tau0
        eetaT[i] = tau0.eeta[i] + exp(alpha0) * (T0[i] - tau0.eeta[i])
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    eetaT[i] = runif(n = 1, min = 0, max = eetaT[i])
  }
  is.observed.eeta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(eetaT[i] > 700){
      eetaT[i] = 700
      is.observed.eeta[i] = FALSE
    }
  }
  
  datlist.eeta = list(N = npatient,
                     time = eetaT,
                     is.observed = is.observed.eeta,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  eetaModel = jags.parallel(data = datlist.eeta,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'),
                           n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                           n.thin = 1, model.file = modelfile.etaModel)
  eeta.sim = cbind(eeta.sim, eetaModel$BUGSoutput[10]$summary[c(1:3, 5:6,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.eeta)){
    if(trueVal.eeta[i] >= eetaModel$BUGSoutput[10]$summary[eeta.order[i], 3]){
      if(trueVal.eeta[i] <= eetaModel$BUGSoutput[10]$summary[eeta.order[i], 7]){
        eeta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(eetaModel$BUGSoutput[10]$summary[1, 3] >= 0){eeta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(eetaModel$BUGSoutput[10]$summary[1, 7] <= 0){eeta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(eetaModel$BUGSoutput[10]$summary[2, 3] >= 0){eeta.significance[2, isim] = TRUE}
  if(eetaModel$BUGSoutput[10]$summary[2, 7] <= 0){eeta.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(eetaModel$BUGSoutput[10]$summary[3, 7] <= 0){eeta.significance[3, isim] = TRUE}
  if(eetaModel$BUGSoutput[10]$summary[3, 3] >= 0){eeta.wrongsignificance[3, isim] = TRUE}
  
  ##eeta1
  if(eetaModel$BUGSoutput[10]$summary[5, 7] <= 0){eeta.significance[4, isim] = TRUE}
  if(eetaModel$BUGSoutput[10]$summary[5, 3] >= 0){eeta.wrongsignificance[4, isim] = TRUE}
  
  ##eeta2
  if(eetaModel$BUGSoutput[10]$summary[6, 3] >= 0){eeta.significance[5, isim] = TRUE}
  if(eetaModel$BUGSoutput[10]$summary[6, 7] <= 0){eeta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(eetaModel$BUGSoutput[10]$summary[9, 3] >= 0){eeta.significance[6, isim] = TRUE}
  if(eetaModel$BUGSoutput[10]$summary[9, 7] <= 0){eeta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated eeta model
trueVal.eeta = c(alpha0, beta10, beta20, eeta10, eeta20, tau0)
(eeta.mean = apply(eeta.sim, FUN = mean, 1))
(eeta.bias = apply(eeta.sim - trueVal.eeta, FUN = mean, 1))

eeta400 = cbind(trueVal.eeta, eeta.mean, eeta.bias)
eeta.ci = rbind(as.numeric(t.test(eeta.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eeta.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eeta.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eeta.sim[4,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eeta.sim[5,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eeta.sim[6,], conf.level = 0.95)$conf.int))
eeta400 = cbind(eeta400, eeta.ci, apply(eeta.cover, FUN = mean, 1), 
               apply(eeta.significance, FUN = mean, 1),
               apply(eeta.wrongsignificance, FUN = mean, 1))
colnames(eeta400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
eeta400






library(runjags)
datlist.eeta$is.observed = as.integer(datlist.eeta$is.observed)
datlist.eeta$trt = as.integer(datlist.eeta$trt)
test = run.jags(modelstring.etaModel, data=datlist.eeta, 
                monitor=c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'), 
                sample=2000, burnin = 4000, adapt = 1000, thin = 1, n.chains = 3, 
                method='parallel')








########################################################################################
########################################################################################
modelstring.delta.inform = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  alpha ~ dnorm(1.506882, 5)
  beta[1] ~ dnorm(2.004770, 5)
  beta[2] ~ dnorm(-1.807582, 5)
   
  delta[1] ~ dnorm(1.5, 2)
  delta[2] ~ dnorm(-1.5, 2)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(2.553775, 10)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.inform = file.path(dir, 'delta.inform.txt')
writeLines(modelstring.delta.inform, con=modelfile.delta.inform)

##why dnorm(1.5, prec = 1/10)?
ex = rnorm(400, mean = 1.5, sd = sqrt(1))
quantile(ex,probs=c(.025,.975))


##simulation for delta model
set.seed(111)
npatient = 400
nsim = 1000
delta.inform.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.inform.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.inform.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.inform.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  delta.inform = jags.parallel(data = datlist.delta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                             n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                             n.thin = 1, model.file = modelfile.delta.inform)
  delta.inform.sim = cbind(delta.inform.sim, delta.inform$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= delta.inform$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= delta.inform$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.inform.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 3] >= 0){delta.inform.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 7] <= 0){delta.inform.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(delta.inform$BUGSoutput[10]$summary[2, 3] >= 0){delta.inform.significance[2, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[2, 7] <= 0){delta.inform.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(delta.inform$BUGSoutput[10]$summary[3, 7] <= 0){delta.inform.significance[3, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[3, 3] >= 0){delta.inform.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(delta.inform$BUGSoutput[10]$summary[4, 3] >= 0){delta.inform.significance[4, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[4, 7] <= 0){delta.inform.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(delta.inform$BUGSoutput[10]$summary[5, 7] <= 0){delta.inform.significance[5, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[5, 3] >= 0){delta.inform.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(delta.inform$BUGSoutput[10]$summary[9, 3] >= 0){delta.inform.significance[6, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[9, 7] <= 0){delta.inform.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta.inform = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.inform.mean = apply(delta.inform.sim, FUN = mean, 1))
(delta.inform.bias = apply(delta.inform.sim - trueVal.delta.inform, FUN = mean, 1))

delta.inform400 = cbind(trueVal.delta.inform, delta.inform.mean, delta.inform.bias)
delta.inform.ci = rbind(as.numeric(t.test(delta.inform.sim[1,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.inform.sim[2,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.inform.sim[3,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.inform.sim[4,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.inform.sim[5,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.inform.sim[6,], conf.level = 0.95)$conf.int))
delta.inform400 = cbind(delta.inform400, delta.inform.ci, apply(delta.inform.cover, FUN = mean, 1),
                 apply(delta.inform.significance, FUN = mean, 1),
                 apply(delta.inform.wrongsignificance, FUN = mean, 1))
colnames(delta.inform400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta.inform400



########################################################################################
########################################################################################
modelstring.eta.inform = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    eta.form[i] <- eta[1] * x1[i] + eta[2] * x2[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(alpha) + eta.form[i] * exp(alpha)))
                              - mu - alpha - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(alpha) + eta.form[i] * exp(alpha)))
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
  
  alpha ~ dnorm(1.506882, 5)
  beta[1] ~ dnorm(2.004770, 5)
  beta[2] ~ dnorm(-1.807582, 5)
  
  eta[1] ~ dnorm(-1, 4)
  eta[2] ~ dnorm(1, 4)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(2.516386, 10)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.inform = file.path(dir, 'eta.inform.txt')
writeLines(modelstring.eta.inform, con=modelfile.eta.inform)


##simulation for eta model
set.seed(111)
npatient = 400
nsim = 100
eta.inform.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; eta10 = -1; eta20 = 1
trueVal.eta = c(alpha0, beta10, beta20, eta10, eta20, tau0)
eta.inform.cover = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.inform.cover) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.order = c(1:3, 5:6, 9)
eta.inform.significance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.inform.significance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.inform.wrongsignificance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.inform.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('eta current run', isim))
  
  ##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  ##set eta1 = 1.5 and eta2 = 0.2
  eta10 = -1; eta20 = 1
  
  ##treatment Z ~ Bin(n, 0.5)
  ##x1 ~ N(0.6, 0.4)
  ##x2 ~ LN(-0.8, 0.6)
  ##epsilon ~ Gumbel(0, 1)
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  
  ##now determine which patient has treatment effect and change their T
  tau0.eta = tau0 + eta10 * x1 + eta20 * x2
  etaT = T0
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0.eta[i]){
        ##patients with treatment and T > tau0
        etaT[i] = tau0.eta[i] + exp(alpha0) * (T0[i] - tau0.eta[i])
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    etaT[i] = runif(n = 1, min = 0, max = etaT[i])
  }
  is.observed.eta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(etaT[i] > 700){
      etaT[i] = 700
      is.observed.eta[i] = FALSE
    }
  }
  
  datlist.eta = list(N = npatient,
                     time = etaT,
                     is.observed = is.observed.eta,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  eta.inform = jags.parallel(data = datlist.eta,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'),
                           n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                           n.thin = 1, model.file = modelfile.eta.inform)
  eta.inform.sim = cbind(eta.inform.sim, eta.inform$BUGSoutput[10]$summary[c(1:3, 5:6,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.eta)){
    if(trueVal.eta[i] >= eta.inform$BUGSoutput[10]$summary[eta.order[i], 3]){
      if(trueVal.eta[i] <= eta.inform$BUGSoutput[10]$summary[eta.order[i], 7]){
        eta.inform.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(eta.inform$BUGSoutput[10]$summary[1, 3] >= 0){eta.inform.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(eta.inform$BUGSoutput[10]$summary[1, 7] <= 0){eta.inform.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(eta.inform$BUGSoutput[10]$summary[2, 3] >= 0){eta.inform.significance[2, isim] = TRUE}
  if(eta.inform$BUGSoutput[10]$summary[2, 7] <= 0){eta.inform.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(eta.inform$BUGSoutput[10]$summary[3, 7] <= 0){eta.inform.significance[3, isim] = TRUE}
  if(eta.inform$BUGSoutput[10]$summary[3, 3] >= 0){eta.inform.wrongsignificance[3, isim] = TRUE}
  
  ##eta1
  if(eta.inform$BUGSoutput[10]$summary[5, 7] <= 0){eta.inform.significance[4, isim] = TRUE}
  if(eta.inform$BUGSoutput[10]$summary[5, 3] >= 0){eta.inform.wrongsignificance[4, isim] = TRUE}
  
  ##eta2
  if(eta.inform$BUGSoutput[10]$summary[6, 3] >= 0){eta.inform.significance[5, isim] = TRUE}
  if(eta.inform$BUGSoutput[10]$summary[6, 7] <= 0){eta.inform.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(eta.inform$BUGSoutput[10]$summary[9, 3] >= 0){eta.inform.significance[6, isim] = TRUE}
  if(eta.inform$BUGSoutput[10]$summary[9, 7] <= 0){eta.inform.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated eta model
trueVal.eta.inform = c(alpha0, beta10, beta20, eta10, eta20, tau0)
(eta.mean = apply(eta.inform.sim, FUN = mean, 1))
(eta.bias = apply(eta.inform.sim - trueVal.eta, FUN = mean, 1))

eta.inform400 = cbind(trueVal.eta, eta.mean, eta.bias)
eta.ci = rbind(as.numeric(t.test(eta.inform.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.inform.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.inform.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.inform.sim[4,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.inform.sim[5,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.inform.sim[6,], conf.level = 0.95)$conf.int))
eta.inform400 = cbind(eta.inform400, eta.ci, apply(eta.inform.cover, FUN = mean, 1), 
               apply(eta.inform.significance, FUN = mean, 1),
               apply(eta.inform.wrongsignificance, FUN = mean, 1))
colnames(eta.inform400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
eta.inform400






########################################################################################
########################################################################################
########################################################################################
########################################################################################
modelstring.delta.weak = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  alpha ~ dnorm(1.506882, 5)
  beta[1] ~ dnorm(2.004770, 5)
  beta[2] ~ dnorm(-1.807582, 5)
   
  delta[1] ~ dnorm(1.5, 1)
  delta[2] ~ dnorm(-1.5, 1)
  #delta[1] ~ dunif(-3, 3)
  #delta[2] ~ dunif(-3, 3)
  

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(2.553775, 10)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta.weak = file.path(dir, 'delta.weak.txt')
writeLines(modelstring.delta.weak, con=modelfile.delta.weak)

##
ex = rnorm(400, mean = 1.5, sd = sqrt(1))
quantile(ex,probs=c(.025,.975))


##simulation for delta model
set.seed(111)
npatient = 400
nsim = 1000
delta.weak.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.weak.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.weak.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.weak.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.weak.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.weak.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.weak.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  delta.weak = jags.parallel(data = datlist.delta,
                               parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                               n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                               n.thin = 1, model.file = modelfile.delta.weak)
  delta.weak.sim = cbind(delta.weak.sim, delta.weak$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= delta.weak$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= delta.weak$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.weak.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(delta.weak$BUGSoutput[10]$summary[1, 3] >= 0){delta.weak.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(delta.weak$BUGSoutput[10]$summary[1, 7] <= 0){delta.weak.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(delta.weak$BUGSoutput[10]$summary[2, 3] >= 0){delta.weak.significance[2, isim] = TRUE}
  if(delta.weak$BUGSoutput[10]$summary[2, 7] <= 0){delta.weak.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(delta.weak$BUGSoutput[10]$summary[3, 7] <= 0){delta.weak.significance[3, isim] = TRUE}
  if(delta.weak$BUGSoutput[10]$summary[3, 3] >= 0){delta.weak.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(delta.weak$BUGSoutput[10]$summary[4, 3] >= 0){delta.weak.significance[4, isim] = TRUE}
  if(delta.weak$BUGSoutput[10]$summary[4, 7] <= 0){delta.weak.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(delta.weak$BUGSoutput[10]$summary[5, 7] <= 0){delta.weak.significance[5, isim] = TRUE}
  if(delta.weak$BUGSoutput[10]$summary[5, 3] >= 0){delta.weak.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(delta.weak$BUGSoutput[10]$summary[9, 3] >= 0){delta.weak.significance[6, isim] = TRUE}
  if(delta.weak$BUGSoutput[10]$summary[9, 7] <= 0){delta.weak.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta.weak = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.weak.mean = apply(delta.weak.sim, FUN = mean, 1))
(delta.weak.bias = apply(delta.weak.sim - trueVal.delta.weak, FUN = mean, 1))

delta.weak400 = cbind(trueVal.delta.weak, delta.weak.mean, delta.weak.bias)
delta.weak.ci = rbind(as.numeric(t.test(delta.weak.sim[1,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.weak.sim[2,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.weak.sim[3,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.weak.sim[4,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.weak.sim[5,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.weak.sim[6,], conf.level = 0.95)$conf.int))
delta.weak400 = cbind(delta.weak400, delta.weak.ci, apply(delta.weak.cover, FUN = mean, 1),
                        apply(delta.weak.significance, FUN = mean, 1),
                        apply(delta.weak.wrongsignificance, FUN = mean, 1))
colnames(delta.weak400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta.weak400


delta.weak$BUGSoutput$summary[,c(1,3,7)]
deltaModel$BUGSoutput$summary[,c(1,3,7)]



########################################################################################
########################################################################################
modelstring.delta2.weak = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  alpha ~ dnorm(1.506882, 5)
  beta[1] ~ dnorm(2.004770, 5)
  beta[2] ~ dnorm(-1.807582, 5)
   
  #delta[1] ~ dnorm(1.5, 1)
  #delta[2] ~ dnorm(-1.5, 1)
  delta[1] ~ dunif(-3, 3)
  delta[2] ~ dunif(-3, 3)
  

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(2.553775, 10)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.delta2.weak = file.path(dir, 'delta2.weak.txt')
writeLines(modelstring.delta2.weak, con=modelfile.delta2.weak)

##
ex = rnorm(400, mean = 1.5, sd = sqrt(1))
quantile(ex,probs=c(.025,.975))


##simulation for delta model
set.seed(111)
npatient = 400
nsim = 1000
delta2.weak.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta2.weak.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta2.weak.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta2.weak.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta2.weak.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta2.weak.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta2.weak.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  delta2.weak = jags.parallel(data = datlist.delta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                             n.chains = 3, n.burnin = 5000, n.iter = 6000, 
                             n.thin = 1, model.file = modelfile.delta2.weak)
  delta2.weak.sim = cbind(delta2.weak.sim, delta2.weak$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= delta2.weak$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= delta2.weak$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta2.weak.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(delta2.weak$BUGSoutput[10]$summary[1, 3] >= 0){delta2.weak.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(delta2.weak$BUGSoutput[10]$summary[1, 7] <= 0){delta2.weak.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(delta2.weak$BUGSoutput[10]$summary[2, 3] >= 0){delta2.weak.significance[2, isim] = TRUE}
  if(delta2.weak$BUGSoutput[10]$summary[2, 7] <= 0){delta2.weak.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(delta2.weak$BUGSoutput[10]$summary[3, 7] <= 0){delta2.weak.significance[3, isim] = TRUE}
  if(delta2.weak$BUGSoutput[10]$summary[3, 3] >= 0){delta2.weak.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(delta2.weak$BUGSoutput[10]$summary[4, 3] >= 0){delta2.weak.significance[4, isim] = TRUE}
  if(delta2.weak$BUGSoutput[10]$summary[4, 7] <= 0){delta2.weak.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(delta2.weak$BUGSoutput[10]$summary[5, 7] <= 0){delta2.weak.significance[5, isim] = TRUE}
  if(delta2.weak$BUGSoutput[10]$summary[5, 3] >= 0){delta2.weak.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(delta2.weak$BUGSoutput[10]$summary[9, 3] >= 0){delta2.weak.significance[6, isim] = TRUE}
  if(delta2.weak$BUGSoutput[10]$summary[9, 7] <= 0){delta2.weak.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta2.weak = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta2.weak.mean = apply(delta2.weak.sim, FUN = mean, 1))
(delta2.weak.bias = apply(delta2.weak.sim - trueVal.delta2.weak, FUN = mean, 1))

delta2.weak400 = cbind(trueVal.delta2.weak, delta2.weak.mean, delta2.weak.bias)
delta2.weak.ci = rbind(as.numeric(t.test(delta2.weak.sim[1,], conf.level = 0.95)$conf.int),
                      as.numeric(t.test(delta2.weak.sim[2,], conf.level = 0.95)$conf.int),
                      as.numeric(t.test(delta2.weak.sim[3,], conf.level = 0.95)$conf.int),
                      as.numeric(t.test(delta2.weak.sim[4,], conf.level = 0.95)$conf.int),
                      as.numeric(t.test(delta2.weak.sim[5,], conf.level = 0.95)$conf.int),
                      as.numeric(t.test(delta2.weak.sim[6,], conf.level = 0.95)$conf.int))
delta2.weak400 = cbind(delta2.weak400, delta2.weak.ci, apply(delta2.weak.cover, FUN = mean, 1),
                      apply(delta2.weak.significance, FUN = mean, 1),
                      apply(delta2.weak.wrongsignificance, FUN = mean, 1))
colnames(delta2.weak400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta2.weak400









########################################################################################
########################################################################################
modelstring.eta.weak = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    eta.form[i] <- eta[1] * x1[i] + eta[2] * x2[i]
    
    ##when treatment = T and Time_i >= tau + eta*Z        
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau - eta.form[i] 
                              + tau * exp(alpha) + eta.form[i] * exp(alpha)))
                              - mu - alpha - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau - eta.form[i]
                    + tau * exp(alpha) + eta.form[i] * exp(alpha)))
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
  
  alpha ~ dnorm(1.506882, 5)
  beta[1] ~ dnorm(2.004770, 5)
  beta[2] ~ dnorm(-1.807582, 5)
  
  #eta[1] ~ dnorm(-1, 2)
  #eta[2] ~ dnorm(1, 2)
  eta[1] ~ dunif(-2.5, 2.5)
  eta[2] ~ dunif(-2.5, 2.5)
  #eta[1] ~ dunif(-0.8*tau, 0.8*tau) 
  #eta[2] ~ dunif(-0.8*tau, 0.8*tau) 

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(2.5, 10)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.eta.weak = file.path(dir, 'eta.weak.txt')
writeLines(modelstring.eta.weak, con=modelfile.eta.weak)

quantile(rnorm(400, mean = 0, sd = sqrt(2)),probs=c(.025,.975))

##simulation for eta model
set.seed(111)
npatient = 400
nsim = 100
eta.weak.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; eta10 = -1; eta20 = 1
trueVal.eta = c(alpha0, beta10, beta20, eta10, eta20, tau0)
eta.weak.cover = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.weak.cover) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.order = c(1:3, 5:6, 9)
eta.weak.significance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.weak.significance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
eta.weak.wrongsignificance = matrix(FALSE, nrow = length(trueVal.eta), ncol = nsim)
rownames(eta.weak.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'eta1', 'eta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('eta current run', isim))
  
  ##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  ##set eta1 = 1.5 and eta2 = 0.2
  eta10 = -1; eta20 = 1
  
  ##treatment Z ~ Bin(n, 0.5)
  ##x1 ~ N(0.6, 0.4)
  ##x2 ~ LN(-0.8, 0.6)
  ##epsilon ~ Gumbel(0, 1)
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  
  ##now determine which patient has treatment effect and change their T
  tau0.eta = tau0 + eta10 * x1 + eta20 * x2
  etaT = T0
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0.eta[i]){
        ##patients with treatment and T > tau0
        etaT[i] = tau0.eta[i] + exp(alpha0) * (T0[i] - tau0.eta[i])
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    etaT[i] = runif(n = 1, min = 0, max = etaT[i])
  }
  is.observed.eta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(etaT[i] > 700){
      etaT[i] = 700
      is.observed.eta[i] = FALSE
    }
  }
  
  datlist.eta = list(N = npatient,
                     time = etaT,
                     is.observed = ifelse(is.observed.eta, 1, 0),
                     trt = ifelse(Z, 1, 0),
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  eta.weak = jags.parallel(data = datlist.eta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'),
                             n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                             n.thin = 1, model.file = modelfile.eta.weak)
  eta.weak.sim = cbind(eta.weak.sim, eta.weak$BUGSoutput[10]$summary[c(1:3, 5:6,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.eta)){
    if(trueVal.eta[i] >= eta.weak$BUGSoutput[10]$summary[eta.order[i], 3]){
      if(trueVal.eta[i] <= eta.weak$BUGSoutput[10]$summary[eta.order[i], 7]){
        eta.weak.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(eta.weak$BUGSoutput[10]$summary[1, 3] >= 0){eta.weak.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(eta.weak$BUGSoutput[10]$summary[1, 7] <= 0){eta.weak.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(eta.weak$BUGSoutput[10]$summary[2, 3] >= 0){eta.weak.significance[2, isim] = TRUE}
  if(eta.weak$BUGSoutput[10]$summary[2, 7] <= 0){eta.weak.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(eta.weak$BUGSoutput[10]$summary[3, 7] <= 0){eta.weak.significance[3, isim] = TRUE}
  if(eta.weak$BUGSoutput[10]$summary[3, 3] >= 0){eta.weak.wrongsignificance[3, isim] = TRUE}
  
  ##eta1
  if(eta.weak$BUGSoutput[10]$summary[5, 7] <= 0){eta.weak.significance[4, isim] = TRUE}
  if(eta.weak$BUGSoutput[10]$summary[5, 3] >= 0){eta.weak.wrongsignificance[4, isim] = TRUE}
  
  ##eta2
  if(eta.weak$BUGSoutput[10]$summary[6, 3] >= 0){eta.weak.significance[5, isim] = TRUE}
  if(eta.weak$BUGSoutput[10]$summary[6, 7] <= 0){eta.weak.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(eta.weak$BUGSoutput[10]$summary[9, 3] >= 0){eta.weak.significance[6, isim] = TRUE}
  if(eta.weak$BUGSoutput[10]$summary[9, 7] <= 0){eta.weak.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated eta model
trueVal.eta.weak = c(alpha0, beta10, beta20, eta10, eta20, tau0)
(eta.mean = apply(eta.weak.sim, FUN = mean, 1))
(eta.bias = apply(eta.weak.sim - trueVal.eta, FUN = mean, 1))

eta.weak400 = cbind(trueVal.eta, eta.mean, eta.bias)
eta.ci = rbind(as.numeric(t.test(eta.weak.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.weak.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.weak.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.weak.sim[4,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.weak.sim[5,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(eta.weak.sim[6,], conf.level = 0.95)$conf.int))
eta.weak400 = cbind(eta.weak400, eta.ci, apply(eta.weak.cover, FUN = mean, 1), 
                      apply(eta.weak.significance, FUN = mean, 1),
                      apply(eta.weak.wrongsignificance, FUN = mean, 1))
colnames(eta.weak400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
eta.weak400


aeta.weak = autorun.jags(model = modelfile.eta.weak, data = datlist.eta, 
                         monitor = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 1, method = 'parallel')
##n=800, 29084 runs, 11 minutes taken
##n=1600, 59 minutes

aeta.weak2 = autorun.jags(model = modelfile.eta.weak, data = datlist.eta, 
                         monitor = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'), 
                         n.chains = 3, startburnin = 5000, startsample = 6000, 
                         thin = 1, method = 'parallel')
##n=400, 2.6 minutes

tau.weak2 = autorun.jags(model = modelfile.tauModel, data = datlist.eta, 
                          monitor = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'), 
                          n.chains = 3, startburnin = 5000, startsample = 6000, 
                          thin = 1, method = 'parallel')

########################################################################################
########################################################################################
########################################################################################
########################################################################################
modelstring.etaNModel = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    eta.form[i] <- eta[1] * x1[i] + eta[2] * x2[i]
    
    ##when treatment = T and Time_i >= tau *exp(eta*Z)      
    epsilon.form.tau[i] <- b * (log(abs(time[i] - tau * exp(eta.form[i])
                              + tau * exp(alpha + eta.form[i])))
                              - mu - alpha - beta.form[i])
    S.tau[i] <- exp(-exp(epsilon.form.tau[i]))
    h.tau[i] <- b / (abs(time[i] - tau * exp(eta.form[i])
                    + tau * exp(alpha + eta.form[i])))
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
  
  alpha ~ dnorm(1.506882, 5)
  beta[1] ~ dnorm(2.004770, 5)
  beta[2] ~ dnorm(-1.807582, 5)
  
  #for(i in 1:2){
  #  eta[i] ~ dnorm(0, 0.1)
  #}
  eta[1] ~ dnorm(0.5, 16)
  eta[2] ~ dnorm(-0.5, 16)

  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0.169534, 0.5)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dnorm(2.516386, 10)I(0,)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.etaNModel = file.path(dir, 'etaNModel.txt')
writeLines(modelstring.etaNModel, con=modelfile.etaNModel)


##simulation for etaN model
set.seed(111)
npatient = 800
nsim = 100
etaN.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; etaN10 = 0.5; etaN20 = -0.5
trueVal.etaN = c(alpha0, beta10, beta20, etaN10, etaN20, tau0)
etaN.cover = matrix(FALSE, nrow = length(trueVal.etaN), ncol = nsim)
rownames(etaN.cover) = c('alpha', 'beta1', 'beta2', 'etaN1', 'etaN2', 'tau')
etaN.order = c(1:3, 5:6, 9)
etaN.significance = matrix(FALSE, nrow = length(trueVal.etaN), ncol = nsim)
rownames(etaN.significance) = c('alpha', 'beta1', 'beta2', 'etaN1', 'etaN2', 'tau')
etaN.wrongsignificance = matrix(FALSE, nrow = length(trueVal.etaN), ncol = nsim)
rownames(etaN.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'etaN1', 'etaN2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('etaN current run', isim))
  
  ##set true alpha = 1.5, tau = 2.5, beta1 = 2, beta2 = 1.8
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  ##set etaN1 = 1.5 and etaN2 = 0.2
  etaN10 = 0.5; etaN20 = -0.5
  
  ##treatment Z ~ Bin(n, 0.5)
  ##x1 ~ N(0.6, 0.4)
  ##x2 ~ LN(-0.8, 0.6)
  ##epsilon ~ Gumbel(0, 1)
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  
  ##now determine which patient has treatment effect and change their T
  tau0.etaN = tau0 * exp(etaN10 * x1 + etaN20 * x2)
  etaNT = T0
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0.etaN[i]){
        ##patients with treatment and T > tau0
        etaNT[i] = tau0.etaN[i] + exp(alpha0) * (T0[i] - tau0.etaN[i])
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    etaNT[i] = runif(n = 1, min = 0, max = etaNT[i])
  }
  is.observed.etaN = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(etaNT[i] > 700){
      etaNT[i] = 700
      is.observed.etaN[i] = FALSE
    }
  }
  
  datlist.etaN = list(N = npatient,
                      time = etaNT,
                      is.observed = is.observed.etaN,
                      trt = Z,
                      x1 = x1,
                      x2 = x2,
                      pi = pi)
  etaNModel = jags.parallel(data = datlist.etaN,
                            parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'eta'),
                            n.chains = 3, n.burnin = 7000, n.iter = 8000, 
                            n.thin = 1, model.file = modelfile.etaNModel)
  etaN.sim = cbind(etaN.sim, etaNModel$BUGSoutput[10]$summary[c(1:3, 5:6,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.etaN)){
    if(trueVal.etaN[i] >= etaNModel$BUGSoutput[10]$summary[etaN.order[i], 3]){
      if(trueVal.etaN[i] <= etaNModel$BUGSoutput[10]$summary[etaN.order[i], 7]){
        etaN.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(etaNModel$BUGSoutput[10]$summary[1, 3] >= 0){etaN.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(etaNModel$BUGSoutput[10]$summary[1, 7] <= 0){etaN.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(etaNModel$BUGSoutput[10]$summary[2, 3] >= 0){etaN.significance[2, isim] = TRUE}
  if(etaNModel$BUGSoutput[10]$summary[2, 7] <= 0){etaN.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(etaNModel$BUGSoutput[10]$summary[3, 7] <= 0){etaN.significance[3, isim] = TRUE}
  if(etaNModel$BUGSoutput[10]$summary[3, 3] >= 0){etaN.wrongsignificance[3, isim] = TRUE}
  
  ##etaN1
  if(etaNModel$BUGSoutput[10]$summary[5, 7] <= 0){etaN.significance[4, isim] = TRUE}
  if(etaNModel$BUGSoutput[10]$summary[5, 3] >= 0){etaN.wrongsignificance[4, isim] = TRUE}
  
  ##etaN2
  if(etaNModel$BUGSoutput[10]$summary[6, 3] >= 0){etaN.significance[5, isim] = TRUE}
  if(etaNModel$BUGSoutput[10]$summary[6, 7] <= 0){etaN.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(etaNModel$BUGSoutput[10]$summary[9, 3] >= 0){etaN.significance[6, isim] = TRUE}
  if(etaNModel$BUGSoutput[10]$summary[9, 7] <= 0){etaN.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated etaN model
trueVal.etaN = c(alpha0, beta10, beta20, etaN10, etaN20, tau0)
(etaN.mean = apply(etaN.sim, FUN = mean, 1))
(etaN.bias = apply(etaN.sim - trueVal.etaN, FUN = mean, 1))

etaN400 = cbind(trueVal.etaN, etaN.mean, etaN.bias)
etaN.ci = rbind(as.numeric(t.test(etaN.sim[1,], conf.level = 0.95)$conf.int),
                as.numeric(t.test(etaN.sim[2,], conf.level = 0.95)$conf.int),
                as.numeric(t.test(etaN.sim[3,], conf.level = 0.95)$conf.int),
                as.numeric(t.test(etaN.sim[4,], conf.level = 0.95)$conf.int),
                as.numeric(t.test(etaN.sim[5,], conf.level = 0.95)$conf.int),
                as.numeric(t.test(etaN.sim[6,], conf.level = 0.95)$conf.int))
etaN400 = cbind(etaN400, etaN.ci, apply(etaN.cover, FUN = mean, 1), 
                apply(etaN.significance, FUN = mean, 1),
                apply(etaN.wrongsignificance, FUN = mean, 1))
colnames(etaN400) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
etaN400




##try this formulation for eta model
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
##set eta1 = 1.5 and eta2 = 0.2
eta10 = 0.5; eta20 = -0.5

##treatment Z ~ Bin(n, 0.5)
##x1 ~ N(0.6, 0.4)
##x2 ~ LN(-0.8, 0.6)
##epsilon ~ Gumbel(0, 1)
Z = rbernoulli(n=npatient, p=0.5)
x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))

tau0.eta = tau0*exp(eta10 * x1 + eta20 * x2)
quantile(tau0.eta, probs=c(.025,.975))










########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

##Simulation study
##n=400; 100 times each

##Tau Model
set.seed(111)
npatient = 400
nsim = 1000
tau.sim = c()
tau.mu = rep(0, times = nsim)
trueVal.tau = c(alpha0, beta10, beta20, tau0)
tau.cover = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.cover) = c('alpha', 'beta1', 'beta2', 'tau')
tau.order = c(1, 2, 3, 7)
tau.significance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.significance) = c('alpha', 'beta1', 'beta2', 'tau')
tau.wrongsignificance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('tau current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  tauT = T0
  
  ##now determine which patient has treatment effect and change their T
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        tauT[i] = tau0 + exp(alpha0) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    tauT[i] = runif(n = 1, min = 0, max = tauT[i])
  }
  is.observed.tau = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(tauT[i] > 700){
      tauT[i] = 700
      is.observed.tau[i] = FALSE
    }
  }
  
  datlist.tau = list(N = npatient,
                     time = tauT,
                     is.observed = is.observed.tau,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  tauModel = jags.parallel(data = datlist.tau,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                           n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                           n.thin = 1, model.file = modelfile.tauModel)
  tau.sim = cbind(tau.sim, tauModel$BUGSoutput[10]$summary[c(1:3,7),1])
  
  ##check if this run covers the true coef
  for(i in 1:4){
    if(trueVal.tau[i] >= tauModel$BUGSoutput[10]$summary[tau.order[i], 3]){
      if(trueVal.tau[i] <= tauModel$BUGSoutput[10]$summary[tau.order[i], 7]){
        tau.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(tauModel$BUGSoutput[10]$summary[1, 3] >= 0){tau.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(tauModel$BUGSoutput[10]$summary[1, 7] <= 0){tau.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(tauModel$BUGSoutput[10]$summary[2, 3] >= 0){tau.significance[2, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[2, 7] <= 0){tau.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(tauModel$BUGSoutput[10]$summary[3, 7] <= 0){tau.significance[3, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[3, 3] >= 0){tau.wrongsignificance[3, isim] = TRUE}
  
  ##tau
  if(tauModel$BUGSoutput[10]$summary[7, 3] >= 0){tau.significance[4, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[7, 7] <= 0){tau.wrongsignificance[4, isim] = TRUE}
  
  ##record mu for setting the prior of delta models
  tau.mu[isim] = tauModel$BUGSoutput[10]$summary[5,1]
}
##mean of the simulation of tau model
(tau.mean = apply(tau.sim, FUN = mean, 1))

##mean of the bias of the simulated tau model
trueVal.tau = c(alpha0, beta10, beta20, tau0)
(tau.bias = apply(tau.sim - trueVal.tau, FUN = mean, 1))

mean(tau.mu)

##result
tau400.1000 = cbind(trueVal.tau, tau.mean, tau.bias)
tau.ci = rbind(as.numeric(t.test(tau.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[4,], conf.level = 0.95)$conf.int))
tau400.1000 = cbind(tau400.1000, tau.ci, apply(tau.cover, FUN = mean, 1),
               apply(tau.significance, FUN = mean, 1),
               apply(tau.wrongsignificance, FUN = mean, 1))
tau.sim400 = tau.sim
colnames(tau400.1000) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
tau400.1000


########################################################################################
########################################################################################
##simulation for delta model
set.seed(111)
npatient = 400
nsim = 1000
delta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  deltaModel = jags.parallel(data = datlist.delta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                             n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                             n.thin = 1, model.file = modelfile.deltaModel)
  delta.sim = cbind(delta.sim, deltaModel$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= deltaModel$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= deltaModel$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 3] >= 0){delta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 7] <= 0){delta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(deltaModel$BUGSoutput[10]$summary[2, 3] >= 0){delta.significance[2, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[2, 7] <= 0){delta.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(deltaModel$BUGSoutput[10]$summary[3, 7] <= 0){delta.significance[3, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[3, 3] >= 0){delta.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(deltaModel$BUGSoutput[10]$summary[4, 3] >= 0){delta.significance[4, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[4, 7] <= 0){delta.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(deltaModel$BUGSoutput[10]$summary[5, 7] <= 0){delta.significance[5, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[5, 3] >= 0){delta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(deltaModel$BUGSoutput[10]$summary[9, 3] >= 0){delta.significance[6, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[9, 7] <= 0){delta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.mean = apply(delta.sim, FUN = mean, 1))
(delta.bias = apply(delta.sim - trueVal.delta, FUN = mean, 1))

delta400.1000 = cbind(trueVal.delta, delta.mean, delta.bias)
delta.ci = rbind(as.numeric(t.test(delta.sim[1,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[2,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[3,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[4,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[5,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[6,], conf.level = 0.95)$conf.int))
delta400.1000 = cbind(delta400.1000, delta.ci, apply(delta.cover, FUN = mean, 1),
                 apply(delta.significance, FUN = mean, 1),
                 apply(delta.wrongsignificance, FUN = mean, 1))
colnames(delta400.1000) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta400.1000







########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

##Simulation study
##n=100; 1000 times each

##Tau Model
set.seed(123)
npatient = 100
nsim = 1000
tau.sim = c()
tau.mu = rep(0, times = nsim)
trueVal.tau = c(alpha0, beta10, beta20, tau0)
tau.cover = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.cover) = c('alpha', 'beta1', 'beta2', 'tau')
tau.order = c(1, 2, 3, 7)
tau.significance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.significance) = c('alpha', 'beta1', 'beta2', 'tau')
tau.wrongsignificance = matrix(FALSE, nrow = length(trueVal.tau), ncol = nsim)
rownames(tau.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('tau current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8
  
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  tauT = T0
  
  ##now determine which patient has treatment effect and change their T
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        tauT[i] = tau0 + exp(alpha0) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    tauT[i] = runif(n = 1, min = 0, max = tauT[i])
  }
  is.observed.tau = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(tauT[i] > 700){
      tauT[i] = 700
      is.observed.tau[i] = FALSE
    }
  }
  
  datlist.tau = list(N = npatient,
                     time = tauT,
                     is.observed = is.observed.tau,
                     trt = Z,
                     x1 = x1,
                     x2 = x2,
                     pi = pi)
  tauModel = jags.parallel(data = datlist.tau,
                           parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha'),
                           n.chains = 3, n.burnin = 3000, n.iter = 5000, 
                           n.thin = 1, model.file = modelfile.tauModel)
  tau.sim = cbind(tau.sim, tauModel$BUGSoutput[10]$summary[c(1:3,7),1])
  
  ##check if this run covers the true coef
  for(i in 1:4){
    if(trueVal.tau[i] >= tauModel$BUGSoutput[10]$summary[tau.order[i], 3]){
      if(trueVal.tau[i] <= tauModel$BUGSoutput[10]$summary[tau.order[i], 7]){
        tau.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(tauModel$BUGSoutput[10]$summary[1, 3] >= 0){tau.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(tauModel$BUGSoutput[10]$summary[1, 7] <= 0){tau.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(tauModel$BUGSoutput[10]$summary[2, 3] >= 0){tau.significance[2, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[2, 7] <= 0){tau.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(tauModel$BUGSoutput[10]$summary[3, 7] <= 0){tau.significance[3, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[3, 3] >= 0){tau.wrongsignificance[3, isim] = TRUE}
  
  ##tau
  if(tauModel$BUGSoutput[10]$summary[7, 3] >= 0){tau.significance[4, isim] = TRUE}
  if(tauModel$BUGSoutput[10]$summary[7, 7] <= 0){tau.wrongsignificance[4, isim] = TRUE}
  
  ##record mu for setting the prior of delta models
  tau.mu[isim] = tauModel$BUGSoutput[10]$summary[5,1]
}
##mean of the simulation of tau model
(tau.mean = apply(tau.sim, FUN = mean, 1))

##mean of the bias of the simulated tau model
trueVal.tau = c(alpha0, beta10, beta20, tau0)
(tau.bias = apply(tau.sim - trueVal.tau, FUN = mean, 1))

mean(tau.mu)

##result
tau100 = cbind(trueVal.tau, tau.mean, tau.bias)
tau.ci = rbind(as.numeric(t.test(tau.sim[1,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[2,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[3,], conf.level = 0.95)$conf.int),
               as.numeric(t.test(tau.sim[4,], conf.level = 0.95)$conf.int))
tau100 = cbind(tau100, tau.ci, apply(tau.cover, FUN = mean, 1),
               apply(tau.significance, FUN = mean, 1),
               apply(tau.wrongsignificance, FUN = mean, 1))
tau.sim100 = tau.sim
colnames(tau100) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
tau100






set.seed(123)
npatient = 100
nsim = 1000
delta.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  deltaModel = jags.parallel(data = datlist.delta,
                             parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                             n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                             n.thin = 1, model.file = modelfile.deltaModel)
  delta.sim = cbind(delta.sim, deltaModel$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= deltaModel$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= deltaModel$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 3] >= 0){delta.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(deltaModel$BUGSoutput[10]$summary[1, 7] <= 0){delta.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(deltaModel$BUGSoutput[10]$summary[2, 3] >= 0){delta.significance[2, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[2, 7] <= 0){delta.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(deltaModel$BUGSoutput[10]$summary[3, 7] <= 0){delta.significance[3, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[3, 3] >= 0){delta.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(deltaModel$BUGSoutput[10]$summary[4, 3] >= 0){delta.significance[4, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[4, 7] <= 0){delta.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(deltaModel$BUGSoutput[10]$summary[5, 7] <= 0){delta.significance[5, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[5, 3] >= 0){delta.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(deltaModel$BUGSoutput[10]$summary[9, 3] >= 0){delta.significance[6, isim] = TRUE}
  if(deltaModel$BUGSoutput[10]$summary[9, 7] <= 0){delta.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.mean = apply(delta.sim, FUN = mean, 1))
(delta.bias = apply(delta.sim - trueVal.delta, FUN = mean, 1))

delta100 = cbind(trueVal.delta, delta.mean, delta.bias)
delta.ci = rbind(as.numeric(t.test(delta.sim[1,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[2,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[3,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[4,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[5,], conf.level = 0.95)$conf.int),
                 as.numeric(t.test(delta.sim[6,], conf.level = 0.95)$conf.int))
delta100 = cbind(delta100, delta.ci, apply(delta.cover, FUN = mean, 1),
                 apply(delta.significance, FUN = mean, 1),
                 apply(delta.wrongsignificance, FUN = mean, 1))
colnames(delta100) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta100






#####################
#####################

modelstring.info1 = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  alpha ~ dnorm(0, 0.001)
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  delta[1] ~ dunif(-3, 3)
  delta[2] ~ dunif(-3, 3)
  
  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5,5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.info1 = file.path(dir, 'delta.info1.txt')
writeLines(modelstring.info1, con=modelfile.info1)


##simulation for delta model
set.seed(111)
npatient = 100
nsim = 1000
delta.inform.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.inform.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.inform.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.inform.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  delta.inform = jags.parallel(data = datlist.delta,
                               parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                               n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                               n.thin = 1, model.file = modelfile.info1)
  delta.inform.sim = cbind(delta.inform.sim, delta.inform$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= delta.inform$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= delta.inform$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.inform.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 3] >= 0){delta.inform.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 7] <= 0){delta.inform.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(delta.inform$BUGSoutput[10]$summary[2, 3] >= 0){delta.inform.significance[2, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[2, 7] <= 0){delta.inform.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(delta.inform$BUGSoutput[10]$summary[3, 7] <= 0){delta.inform.significance[3, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[3, 3] >= 0){delta.inform.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(delta.inform$BUGSoutput[10]$summary[4, 3] >= 0){delta.inform.significance[4, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[4, 7] <= 0){delta.inform.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(delta.inform$BUGSoutput[10]$summary[5, 7] <= 0){delta.inform.significance[5, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[5, 3] >= 0){delta.inform.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(delta.inform$BUGSoutput[10]$summary[9, 3] >= 0){delta.inform.significance[6, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[9, 7] <= 0){delta.inform.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta.inform = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.inform.mean = apply(delta.inform.sim, FUN = mean, 1))
(delta.inform.bias = apply(delta.inform.sim - trueVal.delta.inform, FUN = mean, 1))

delta.inform400 = cbind(trueVal.delta.inform, delta.inform.mean, delta.inform.bias)
delta.inform.ci = rbind(as.numeric(t.test(delta.inform.sim[1,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[2,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[3,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[4,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[5,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[6,], conf.level = 0.95)$conf.int))
delta.info1 = cbind(delta.inform400, delta.inform.ci, apply(delta.inform.cover, FUN = mean, 1),
                        apply(delta.inform.significance, FUN = mean, 1),
                        apply(delta.inform.wrongsignificance, FUN = mean, 1))
colnames(delta.info1) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta.info1




################
################

modelstring.info2 = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  alpha ~ dnorm(0, 0.001)
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  delta[1] ~ dnorm(1.5, 1)
  delta[2] ~ dnorm(-1.5, 1)
  
  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5,5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.info2 = file.path(dir, 'delta.info2.txt')
writeLines(modelstring.info2, con=modelfile.info2)


##simulation for delta model
set.seed(111)
npatient = 100
nsim = 1000
delta.inform.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.inform.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.inform.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.inform.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  delta.inform = jags.parallel(data = datlist.delta,
                               parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                               n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                               n.thin = 1, model.file = modelfile.info2)
  delta.inform.sim = cbind(delta.inform.sim, delta.inform$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= delta.inform$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= delta.inform$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.inform.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 3] >= 0){delta.inform.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 7] <= 0){delta.inform.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(delta.inform$BUGSoutput[10]$summary[2, 3] >= 0){delta.inform.significance[2, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[2, 7] <= 0){delta.inform.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(delta.inform$BUGSoutput[10]$summary[3, 7] <= 0){delta.inform.significance[3, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[3, 3] >= 0){delta.inform.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(delta.inform$BUGSoutput[10]$summary[4, 3] >= 0){delta.inform.significance[4, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[4, 7] <= 0){delta.inform.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(delta.inform$BUGSoutput[10]$summary[5, 7] <= 0){delta.inform.significance[5, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[5, 3] >= 0){delta.inform.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(delta.inform$BUGSoutput[10]$summary[9, 3] >= 0){delta.inform.significance[6, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[9, 7] <= 0){delta.inform.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta.inform = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.inform.mean = apply(delta.inform.sim, FUN = mean, 1))
(delta.inform.bias = apply(delta.inform.sim - trueVal.delta.inform, FUN = mean, 1))

delta.inform400 = cbind(trueVal.delta.inform, delta.inform.mean, delta.inform.bias)
delta.inform.ci = rbind(as.numeric(t.test(delta.inform.sim[1,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[2,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[3,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[4,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[5,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[6,], conf.level = 0.95)$conf.int))
delta.info2 = cbind(delta.inform400, delta.inform.ci, apply(delta.inform.cover, FUN = mean, 1),
                        apply(delta.inform.significance, FUN = mean, 1),
                        apply(delta.inform.wrongsignificance, FUN = mean, 1))
colnames(delta.info2) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta.info2




###############
###############

modelstring.info3 = '
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
    beta.form[i] <- beta[1] * x1[i] + beta[2] * x2[i]
    
    delta.form[i] <- delta[1] * x1[i] + delta[2] * x2[i]
    
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
  
  alpha ~ dnorm(0, 0.001)
  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.001)
  }
  delta[1] ~ dnorm(1.5, 0.5)
  delta[2] ~ dnorm(-1.5, 0.5)
  
  ##prior for mu and sigma of the weibull distribution
  mu ~ dnorm(0, 0.001)
  b ~ dgamma(0.001, 0.001)
  sigma <- pow(b, -1)
  tau ~ dunif(0.5,5)
}
'
dir = '/Users/xingz/Library/CloudStorage/OneDrive-Personal/Desktop/Every Old Things/Osaka University/Study/Oak Clinical/R'
modelfile.info3 = file.path(dir, 'delta.info3.txt')
writeLines(modelstring.info3, con=modelfile.info3)


##simulation for delta model
set.seed(111)
npatient = 100
nsim = 1000
delta.inform.sim = c()
alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
trueVal.delta = c(alpha0, beta10, beta20, delta10, delta20, tau0)
delta.inform.cover = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.cover) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.order = c(1:5, 9)
delta.inform.significance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.significance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
delta.inform.wrongsignificance = matrix(FALSE, nrow = length(trueVal.delta), ncol = nsim)
rownames(delta.inform.wrongsignificance) = c('alpha', 'beta1', 'beta2', 'delta1', 'delta2', 'tau')
for(isim in 1:nsim){
  ##announce the current number of run
  print(paste('delta current run', isim))
  
  alpha0 = 1.5; tau0 = 2.5; beta10 = 2; beta20 = -1.8; delta10 = 1.5; delta20 = -1.5
  Z = rbernoulli(n=npatient, p=0.5)
  x1 = rnorm(n = npatient, mean = 0.6, sd = sqrt(0.4))
  x2 = rlnorm(n = npatient, meanlog = -0.8, sdlog = sqrt(0.6))
  epsilon = log(rexp(npatient))
  
  ##calculate T0 = beta10 * x1 + beta20 * x2 + epsilon
  T0 = exp(beta10 * x1 + beta20 * x2 + epsilon)
  deltaT = T0
  
  ##now determine which patient has treatment effect and change their T
  alpha0.delta = alpha0 + delta10 * x1 + delta20 * x2
  for(i in 1:npatient){
    if(Z[i]){
      if(T0[i] > tau0){
        ##patients with treatment and T > tau0
        deltaT[i] = tau0 + exp(alpha0.delta[i]) * (T0[i] - tau0)
      }
    }
  }
  
  ##apply censoring. 25% of the patients will be censored
  for(i in 1:(npatient/4)){
    deltaT[i] = runif(n = 1, min = 0, max = deltaT[i])
  }
  is.observed.delta = c(rep(FALSE, times = npatient/4), rep(TRUE, times = 3*npatient/4))
  
  ##set a trial length maximum at 700
  for(i in 1:npatient){
    if(deltaT[i] > 700){
      deltaT[i] = 700
      is.observed.delta[i] = FALSE
    }
  }
  
  datlist.delta = list(N = npatient,
                       time = deltaT,
                       is.observed = is.observed.delta,
                       trt = Z,
                       x1 = x1,
                       x2 = x2,
                       pi = pi)
  delta.inform = jags.parallel(data = datlist.delta,
                               parameters.to.save = c('beta', 'mu', 'sigma', 'tau', 'alpha', 'delta'),
                               n.chains = 3, n.burnin = 8000, n.iter = 9000, 
                               n.thin = 1, model.file = modelfile.info3)
  delta.inform.sim = cbind(delta.inform.sim, delta.inform$BUGSoutput[10]$summary[c(1:5,9),1])
  
  ##check if this run covers the true coef
  for(i in 1:length(trueVal.delta)){
    if(trueVal.delta[i] >= delta.inform$BUGSoutput[10]$summary[delta.order[i], 3]){
      if(trueVal.delta[i] <= delta.inform$BUGSoutput[10]$summary[delta.order[i], 7]){
        delta.inform.cover[i, isim] = TRUE
      }
    }
  }
  
  ##check whether the estimate is significant
  
  ##alpha lower bound is greater than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 3] >= 0){delta.inform.significance[1, isim] = TRUE}
  ##alpha upper bound is wrongly to be smaller than 0
  if(delta.inform$BUGSoutput[10]$summary[1, 7] <= 0){delta.inform.wrongsignificance[1, isim] = TRUE}
  
  ##beta1
  if(delta.inform$BUGSoutput[10]$summary[2, 3] >= 0){delta.inform.significance[2, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[2, 7] <= 0){delta.inform.wrongsignificance[2, isim] = TRUE}
  
  ##beta2
  if(delta.inform$BUGSoutput[10]$summary[3, 7] <= 0){delta.inform.significance[3, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[3, 3] >= 0){delta.inform.wrongsignificance[3, isim] = TRUE}
  
  ##delta1
  if(delta.inform$BUGSoutput[10]$summary[4, 3] >= 0){delta.inform.significance[4, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[4, 7] <= 0){delta.inform.wrongsignificance[4, isim] = TRUE}
  
  ##delta2
  if(delta.inform$BUGSoutput[10]$summary[5, 7] <= 0){delta.inform.significance[5, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[5, 3] >= 0){delta.inform.wrongsignificance[5, isim] = TRUE}
  
  ##tau
  if(delta.inform$BUGSoutput[10]$summary[9, 3] >= 0){delta.inform.significance[6, isim] = TRUE}
  if(delta.inform$BUGSoutput[10]$summary[9, 7] <= 0){delta.inform.wrongsignificance[6, isim] = TRUE}
}


##mean and bias of the simulated delta model
trueVal.delta.inform = c(alpha0, beta10, beta20, delta10, delta20, tau0)
(delta.inform.mean = apply(delta.inform.sim, FUN = mean, 1))
(delta.inform.bias = apply(delta.inform.sim - trueVal.delta.inform, FUN = mean, 1))

delta.inform400 = cbind(trueVal.delta.inform, delta.inform.mean, delta.inform.bias)
delta.inform.ci = rbind(as.numeric(t.test(delta.inform.sim[1,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[2,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[3,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[4,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[5,], conf.level = 0.95)$conf.int),
                        as.numeric(t.test(delta.inform.sim[6,], conf.level = 0.95)$conf.int))
delta.info3 = cbind(delta.inform400, delta.inform.ci, apply(delta.inform.cover, FUN = mean, 1),
                        apply(delta.inform.significance, FUN = mean, 1),
                        apply(delta.inform.wrongsignificance, FUN = mean, 1))
colnames(delta.info3) = c('True', 'Estimated', 'Difference', '2.5%', '97.5%', '95%CI coverage', 'sig rate', 'false sig')
delta.info3









