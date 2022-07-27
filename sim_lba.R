library(tidyverse)
library(rstan)
library(here)
library(tmvtnorm)
library(psych)
library(rtdists)
library(progress)
library(abind)
library(msm)
library(truncnorm)
library(data.table)
library(bayesplot)

set.seed(666)

source("M3_functions.R")
source("lba.R")

##### TESTING 
test <- matrix(NaN,ncol=4, nrow=10000)

colnames(test) <- c("t0","b","A","k")

for (j in 1:10000){
  test[j,"t0"] = rtnorm(1,lower=0,mean = 2,sd = 1)
  test[j,"b"] = rtnorm(1,lower = 5, upper= Inf, mean=50,sd=10)
  test[j,"A"] = rtnorm(1,lower=0, upper=(test[j,"b"]-10), mean=25, sd=5)
  test[j,"k"] = test[j,"b"]-test[j,"A"]}

sd(test[,"t0"])

hist(test[,"b"])
#####

# Function for Simulating Trialwise Activations Complex Span Model Drift as Link to Activations ----
simActs_CSpan <-simData_CSpan <- function(parmsMMM,respOpts,nRetrievals,SetSize,VP_ID,A,b){
  # extract parms
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    filt <- parmsMMM["f"]
    baseA <- parmsMMM["baseA"]
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    filt <- parmsMMM[,"f"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  # compute acts for response categories
  A_IIP <- conA + genA + baseA
  A_IOP <- genA + baseA
  A_DIP <- filt*(conA + genA) + baseA
  A_DOP <- filt*genA + baseA
  A_NPL <- baseA
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }
  
  # Compute Probability with acts as drift.
  LBA_pars <- LBA_v2(acts, as.vector(respOpts),a=A,b=b)
  colnames(LBA_pars[[1]]) <- c("P_IIP","P_IOP","P_DIP","P_DOP","P_NPL")
  
  # Extract information 
  
  v <- LBA_pars[[2]]
  #v<-do.call(rbind, replicate(nrow(acts), LBA_pars[[2]], simplify=FALSE))
  colnames(v)<- c("v_1","v_2","v_3","v_4","v_5")
  
  # Simulate Choices from calculated Probabilities according to
  # LBA Choice Rule 
  
  simdata <- matrix(NA,ncol = ncol(LBA_pars[[1]]),nrow = nrow(LBA_pars[[1]]))
  
  for(id in 1:nrow(LBA_pars[[1]])){
    simdata[id,] <- t(stats::rmultinom(1,nRetrievals,LBA_pars[[1]][id,]))
  }
  
  colnames(simdata) <- c("C_IIP","C_IOP","C_DIP","C_DOP","C_NPL")
  # 
  # TrialNum <- rep(1:nTrials,each=SetSize)
  # ID <- rep(1:n_subjects, each=nTrials*SetSize)
  # 
  # trials <- cbind(ID,TrialNum,acts,drifts)
  
  VP_ID <- rep(1,each=n_obs)
  Trialnum <- rep(1:nTrials,each=SetSize)
  lba_choices <- cbind(VP_ID,Trialnum,acts,LBA_pars[[1]],v,simdata)
  
  
  
  return(lba_choices)
}

# Simulate Activations on Trial Level for Simulating LBA Data #

n_subjects <- 5
respOpt<-respOpt_Cspan(4,8)
nTrials <- 30
SetSize <- 5
n_choice <- 5
n_obs <- nTrials*SetSize

# Set Range for Parameter Means
range_muC <- c(1,30)
range_muA <- c(0.75,0.9)
range_muF <- c(0,1) # fix to 0.5
t0_mu <- 0.5

eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56


sigC <- c(0.1,0.2)
sigA <- c(0.1,0.2)
sigF <- c(0.1,0.2)
sigB <- c(0.0001, 0.1)

# Set Up Progress Bar

pb <- progress_bar$new(format = "  Simulating Data [:bar] :percent eta: :eta",
                       total = n_subjects, clear = FALSE, width= 80)

# Create Object
LBA_dat <- list()
RT_data <- array(dim=c(n_subjects,1,2,n_obs))

for (i in 1:n_subjects)
{
  
  # Generate Random Parameter Set from randomly drawn hyper_pars ----
  
  # Sample Hyper-Parameter Means with C as fixpoint ----
  relCA <- runif(1, min = range_muA[1],max = range_muA[2])
  Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
  Mean_Apar <- Mean_Cpar*(relCA)
  Mean_Fpar <- runif(1, min =range_muF[1], max= range_muF[2])
  Mean_bpar <- 0.1
  log_mu_f <- log(Mean_Fpar/(1-Mean_Fpar))
  
  hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_f, Mean_bpar)
  
  
  # Sample Variances and Set Covariances----
  
  sig_c <- runif(1, min = sigC[1], max = sigC[2])*Mean_Cpar
  sig_a <- runif(1, min = sigA[1], max = sigA[2])*Mean_Apar 
  sig_f <- runif(1, min = sigF[1], max = sigF[2])
  sig_b <- 0.001
  
  sigs <-c(sig_c,sig_a,sig_f,sig_b)
  Sig <- diag(length(hyper_mus))
  
  Sig[1,1] <- (sig_c)^2
  Sig[2,2] <- (sig_a)^2
  Sig[3,3] <- (sig_f)^2
  Sig[4,4] <- (sig_b)^2
  
  
  # Set Correlations for Parameters ----
  
  # Sampe Covariance Matrix Sigma
  
  omega <- rlkjcorr(1,length(hyper_mus),eta)
  
  # Little Hack for fixing coveriance of b to zer0
  
  omega[4,1:3] = omega[1:3,4] = 0
  
  Sigma <- cor2cov(omega,sigs)
  
  
  # Sample Parameters from MVN ----
  
  theta <- as.matrix(tmvtnorm::rtmvnorm(n=n_obs, mean= hyper_mus, sigma=Sigma,
                                        lower=c(0,0,-Inf,0),upper = c(Inf,Inf,Inf,Inf)))
  
  # Merge Parameters to one Matrix
  colnames(theta) <- c("conA","genA","f","baseA")
  theta[,4] <- 0.1
  theta[,3] <- 1 / (1+exp(-theta[,3]))
  
  
  # Sample t0 for RT Simulation
 
  t0 = rtnorm(1,lower=0,mean = 0.5,sd = 0.5)
  b = rtnorm(1,lower = 0, upper= Inf, mean=0.5,sd=1)
  A = rtnorm(1,lower=0, upper=b, mean=0.5, sd=1)
  k = b-A
  
  
  # Simulate Trial-Level Activations and transform them to driftrates and probabilities using LBA CDF
  
  trial_lvlacts <- simActs_CSpan(theta,as.vector(respOpt),1,SetSize,A = A,b=b)
  
  # Scaling the Driftrates 
  
  v1_scaled <-   trial_lvlacts[,"v_1"] *10
  v2_scaled <-   trial_lvlacts[,"v_2"] *10
  v3_scaled <-   trial_lvlacts[,"v_3"] *10
  v4_scaled <-   trial_lvlacts[,"v_4"] *10
  v5_scaled <-   trial_lvlacts[,"v_5"] *10

  # v1_scaled <-   qtnorm(trial_lvlacts[,"v_1"],mean = 2,sd=1,lower=0)
  # v2_scaled <-   qtnorm(trial_lvlacts[,"v_2"],mean = 2,sd=1,lower=0)
  # v3_scaled <-   qtnorm(trial_lvlacts[,"v_3"],mean = 2,sd=1,lower=0)
  # v4_scaled <-   qtnorm(trial_lvlacts[,"v_4"],mean = 2,sd=1,lower=0)
  # v5_scaled <-   qtnorm(trial_lvlacts[,"v_5"],mean = 2,sd=1,lower=0)
  
  trial_lvlacts <- cbind(trial_lvlacts, v1_scaled,v2_scaled,v3_scaled,v4_scaled,v5_scaled)
  
  # Simulating RT Data 
  
  lba_rt<- matrix(NaN,nrow = n_obs,ncol=2)
  
 
  pred <- as.vector(rtdists::rLBA(n_obs,A=A,b=b,t0=t0, 
                                    mean_v=c(mean(trial_lvlacts[,"v1_scaled"]),
                                             mean(trial_lvlacts[,"v2_scaled"]),
                                             mean(trial_lvlacts[,"v3_scaled"]),
                                             mean(trial_lvlacts[,"v4_scaled"]),
                                             mean(trial_lvlacts[,"v5_scaled"])),
                                    sd_v=1,silent = T,distribution = "norm"))
    
    lba_rt[,1] <- pred$rt
    lba_rt[,2] <- pred$response
    

  full <- cbind(trial_lvlacts,lba_rt)
  
  LBA_dat[[i]] <- list(full, lba_rt)
  
  pb$tick()
  
}


# Merge Subject Arrays for STAN Data ----

RT_data <- array(NaN,dim=c(n_subjects,2,n_obs))

for (i in 1:n_subjects) 
{
  
  RT_data[i,1,] <- LBA_dat[[i]][[2]][,1]
  RT_data[i,2,] <- LBA_dat[[i]][[2]][,2]
  
  
}

LBA_dat[[2]]
# Check Range of RT Distribution, Check frequencies of Response Categories

for (i in 1:n_subjects){
  hist(RT_data[i,1,])
  hist(RT_data[i,2,])
}

# Setting up Stan Data 

init_lba <-function() 
{
  
  list(list(v=cbind(runif(stan.dat$NUM_SUBJ,1,5))))
  list(v_mu = cbind(runif(stan.dat$NUM_SUBJ,5,10)))}



stan.dat <- list(RT=RT_data,
                 TEST_LENGTH=n_obs,
                 NUM_SUBJ=n_subjects,
                 NUM_CHOICES=5)


fit_lba_drift <- stan("LBA/LBA_hier.stan",
                      data = stan.dat,
                      iter=500,chains=4,
                      warmup = 150,
                      cores=4,refresh=100,
                      control=list(max_treedepth=15))

print(fit_lba_drift)

post_drift <- extract(fit_lba_drift,pars="v")
post_pred <- extract(fit_lba_drift,pars="pred")
v_subj <- colMeans(post_drift$v)

v_vp1 <- LBA_dat[[1]][[1]]
v_vp1 <- v_vp1[,23:27]
v_vp1 <- colMeans(v_vp1)

v_vp2 <- LBA_dat[[2]][[1]]
v_vp2 <- v_vp2[,23:27]
v_vp2 <- colMeans(v_vp2)

v_vp3 <- LBA_dat[[3]][[1]]
v_vp3 <- v_vp3[,23:27]
v_vp3 <- colMeans(v_vp3)

v_vp4 <- LBA_dat[[4]][[1]]
v_vp4 <- v_vp4[,23:27]
v_vp4 <- colMeans(v_vp4)

v_vp5 <- LBA_dat[[5]][[1]]
v_vp5 <- v_vp5[,23:27]
v_vp5 <- colMeans(v_vp5)

v_real <- rbind(v_vp1,v_vp2,v_vp3,v_vp4,v_vp5)
colnames(v_real) <- c("v_IIP_zsim","v_IOP_zsim","v_DIP_zsim","v_DIOP_zsim","v_NPL_zsim")
v_est <- v_subj
colnames(v_subj) <- c("v_IIP_zest","v_IOP_zest","v_DIP_zest","v_DIOP_zest","v_NPL_zest")

v_LBA <- as.data.frame(cbind(v_real,v_est))


v_LBA %>% pivot_wider(everything())
  ggplot(., aes(y=Estimate,x=Category, group=Origin,color=Origin)) + geom_point()


ggplot(v_LBA,aes(x=v_LBA[1,1:5],y=v_LBA[1,6:10])) + geom_point()

cor(v_subj[1,],v_real[1,])
cor(v_subj[2,],v_real[2,])
cor(v_subj[3,],v_real[3,])
cor(v_subj[4,],v_real[4,])

traceplot(fit_lba_drift,pars="v")
dim(v_vp1)
traceplot(fit_lba_drift,"v")
saveRDS(fit_lba_drift,"fit_lba_drift_log.RDS")
saveRDS(LBA_dat,"LBA_dat_driftlink.RDS")

#############################SIMULATE DATA For Link to Starting Point Bias of each Category#################################################


bayesplot::mcmc_acf(fit_lba_drift, regex_pars = "v_mu")

# Function for Simulating Trialwise Activations Complex Span Model Drift as Link to Activations ----
simActs_CSpan <-simData_CSpan <- function(parmsMMM,respOpts,nRetrievals,SetSize,VP_ID,A,b){
  # extract parms
  if(is.vector(parmsMMM)){
    conA <- parmsMMM["conA"]
    genA <- parmsMMM["genA"]
    filt <- parmsMMM["f"]
    baseA <- parmsMMM["baseA"]
  }else{
    conA <- parmsMMM[,"conA"]
    genA <- parmsMMM[,"genA"]
    filt <- parmsMMM[,"f"]
    baseA <- parmsMMM[,"baseA"]
  }
  
  # compute acts for response categories
  A_IIP <- conA + genA + baseA
  A_IOP <- genA + baseA
  A_DIP <- filt*(conA + genA) + baseA
  A_DOP <- filt*genA + baseA
  A_NPL <- baseA
  
  # summarize activations
  if(length(A_IIP) != 1){
    acts <- cbind(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }else{
    acts <- c(A_IIP,A_IOP,A_DIP,A_DOP,A_NPL)
  }
  
  # Compute Probability with acts as drift.
  LBA_pars <- LBA_v2(acts, as.vector(respOpts),a=A,b=b)
  colnames(LBA_pars[[1]]) <- c("P_IIP","P_IOP","P_DIP","P_DOP","P_NPL")
  
  # Extract information 
  
  v <- LBA_pars[[2]]
  #v<-do.call(rbind, replicate(nrow(acts), LBA_pars[[2]], simplify=FALSE))
  colnames(v)<- c("v_1","v_2","v_3","v_4","v_5")
  
  # Simulate Choices from calculated Probabilities according to
  # LBA Choice Rule 
  
  simdata <- matrix(NA,ncol = ncol(LBA_pars[[1]]),nrow = nrow(LBA_pars[[1]]))
  
  for(id in 1:nrow(LBA_pars[[1]])){
    simdata[id,] <- t(stats::rmultinom(1,nRetrievals,LBA_pars[[1]][id,]))
  }
  
  colnames(simdata) <- c("C_IIP","C_IOP","C_DIP","C_DOP","C_NPL")
  # 
  # TrialNum <- rep(1:nTrials,each=SetSize)
  # ID <- rep(1:n_subjects, each=nTrials*SetSize)
  # 
  # trials <- cbind(ID,TrialNum,acts,drifts)
  
  VP_ID <- rep(1,each=n_obs)
  Trialnum <- rep(1:nTrials,each=SetSize)
  lba_choices <- cbind(VP_ID,Trialnum,acts,LBA_pars[[1]],v,simdata)
  
  
  
  return(lba_choices)
}

# Simulate Activations on Trial Level for Simulating LBA Data #

n_subjects <- 5
respOpt<-respOpt_Cspan(4,8)
nTrials <- 20
SetSize <- 5
n_choice <- 5
n_obs <- nTrials*SetSize

# Set Range for Parameter Means
range_muC <- c(1,30)
range_muA <- c(0.75,0.9)
range_muF <- c(0,1) # fix to 0.5
t0_mu <- 0.5

eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56


sigC <- c(0.1,0.2)
sigA <- c(0.1,0.2)
sigF <- c(0.1,0.2)
sigB <- c(0.0001, 0.1)

# Set Up Progress Bar

pb <- progress_bar$new(format = "  Simulating Data [:bar] :percent eta: :eta",
                       total = n_subjects, clear = FALSE, width= 80)

# Create Object
LBA_dat <- list()
RT_data <- array(dim=c(n_subjects,1,2,n_obs))

for (i in 1:n_subjects)
{
  
  # Generate Random Parameter Set from randomly drawn hyper_pars ----
  
  # Sample Hyper-Parameter Means with C as fixpoint ----
  relCA <- runif(1, min = range_muA[1],max = range_muA[2])
  Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
  Mean_Apar <- Mean_Cpar*(relCA)
  Mean_Fpar <- runif(1, min =range_muF[1], max= range_muF[2])
  Mean_bpar <- 0.1
  log_mu_f <- log(Mean_Fpar/(1-Mean_Fpar))
  
  hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_f, Mean_bpar)
  
  
  # Sample Variances and Set Covariances----
  
  sig_c <- runif(1, min = sigC[1], max = sigC[2])*Mean_Cpar
  sig_a <- runif(1, min = sigA[1], max = sigA[2])*Mean_Apar 
  sig_f <- runif(1, min = sigF[1], max = sigF[2])
  sig_b <- 0.001
  
  sigs <-c(sig_c,sig_a,sig_f,sig_b)
  Sig <- diag(length(hyper_mus))
  
  Sig[1,1] <- (sig_c)^2
  Sig[2,2] <- (sig_a)^2
  Sig[3,3] <- (sig_f)^2
  Sig[4,4] <- (sig_b)^2
  
  
  # Set Correlations for Parameters ----
  
  # Sampe Covariance Matrix Sigma
  
  omega <- rlkjcorr(1,length(hyper_mus),eta)
  
  # Little Hack for fixing coveriance of b to zer0
  
  omega[4,1:3] = omega[1:3,4] = 0
  
  Sigma <- cor2cov(omega,sigs)
  
  
  # Sample Parameters from MVN ----
  
  theta <- as.matrix(tmvtnorm::rtmvnorm(n=n_obs, mean= hyper_mus, sigma=Sigma,
                                        lower=c(0,0,-Inf,0),upper = c(Inf,Inf,Inf,Inf)))
  
  # Merge Parameters to one Matrix
  colnames(theta) <- c("conA","genA","f","baseA")
  theta[,4] <- 0.1
  theta[,3] <- 1 / (1+exp(-theta[,3]))
  
  
  # Sample t0 for RT Simulation
  
  t0 = rtnorm(1,lower=0,mean = 0.5,sd = 0.5)
  b = rtnorm(1,lower = 0, upper= Inf, mean=0.5,sd=1)
  A = rtnorm(1,lower=0, upper=A, mean=0.5, sd=1)
  k = b-A
  
  
  # Simulate Trial-Level Activations and transform them to driftrates and probabilities using LBA CDF
  
  trial_lvlacts <- simActs_CSpan(theta,as.vector(respOpt),1,SetSize,A = A,b=b)
  
  # Scaling the Driftrates 
  
  v1_scaled <-   trial_lvlacts[,"v_1"] *10
  v2_scaled <-   trial_lvlacts[,"v_2"] *10
  v3_scaled <-   trial_lvlacts[,"v_3"] *10
  v4_scaled <-   trial_lvlacts[,"v_4"] *10
  v5_scaled <-   trial_lvlacts[,"v_5"] *10
  
  # v1_scaled <-   qtnorm(trial_lvlacts[,"v_1"],mean = 2,sd=1,lower=0)
  # v2_scaled <-   qtnorm(trial_lvlacts[,"v_2"],mean = 2,sd=1,lower=0)
  # v3_scaled <-   qtnorm(trial_lvlacts[,"v_3"],mean = 2,sd=1,lower=0)
  # v4_scaled <-   qtnorm(trial_lvlacts[,"v_4"],mean = 2,sd=1,lower=0)
  # v5_scaled <-   qtnorm(trial_lvlacts[,"v_5"],mean = 2,sd=1,lower=0)
  # 
  trial_lvlacts <- cbind(trial_lvlacts, v1_scaled,v2_scaled,v3_scaled,v4_scaled,v5_scaled)
  
  # Simulating RT Data 
  
  lba_rt<- matrix(NaN,nrow = n_obs,ncol=2)
  
  
  pred <- as.vector(rtdists::rLBA(n_obs,A=A,b=b,t0=t0, 
                                  mean_v=c(mean(trial_lvlacts[,"v1_scaled"]),
                                           mean(trial_lvlacts[,"v2_scaled"]),
                                           mean(trial_lvlacts[,"v3_scaled"]),
                                           mean(trial_lvlacts[,"v4_scaled"]),
                                           mean(trial_lvlacts[,"v5_scaled"])),
                                  sd_v=1,silent = T,distribution = "norm"))
  
  lba_rt[,1] <- pred$rt
  lba_rt[,2] <- pred$response
  
  
  full <- cbind(trial_lvlacts,lba_rt)
  
  LBA_dat[[i]] <- list(full, lba_rt)
  
  pb$tick()
  
}









