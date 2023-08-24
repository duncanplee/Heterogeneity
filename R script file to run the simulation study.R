#################################################
#### Run the different simulation study scenarios 
#################################################
#### Set the seed to make the results reproducible
set.seed(1)


#################################
#### Load libraries and functions
#################################
library(meta)
library(metafor)


source('R simulation function.R')



######################################
#### Specify the simulation quantities
######################################
#### Number of simulations per scenario
N.sim <- 1000


#### Number of studies in each region
N.A <- c(3, 5, 10, 20, 30)
N.B <- c(3, 5, 10, 20, 30)
n.A <- length(N.A)
n.B <- length(N.B)


#### True relative risks
RR.A <- 1.005
RR.B <- 1.005  # Possible values are 1.005, 1.009, 1.013


#### True estimated standard errors
SE.A <- 0.005 
SE.B <- 0.005 
  

####  SD of the beta estimates capturing within region heterogeneity
T.A <- 0.005  # Possible values are 0.0025, 0.005 and 0.01
T.B <- 0.005  # Possible values are 0.0025, 0.005 and 0.01


####  SD of the standard error estimates
V.A <- SE.A / 3  
V.B <- SE.B / 3  



##################
#### Run the study
##################
#### Create a results matrix to store the results across all 25 combinations of N.A and N.B
n.temp <- n.A * n.B
simres <- data.frame(N.A=rep(NA, n.temp), N.B=rep(NA, n.temp), 
                      bias.global=rep(NA, n.temp), bias.SG=rep(NA, n.temp), bias.MR=rep(NA, n.temp),
                      rmse.global=rep(NA, n.temp), rmse.SG=rep(NA, n.temp), rmse.MR=rep(NA, n.temp),
                      coverage.global=rep(NA, n.temp), coverage.SG=rep(NA, n.temp), coverage.MR=rep(NA, n.temp),
                      width.global=rep(NA, n.temp), width.SG=rep(NA, n.temp), width.MR=rep(NA, n.temp), 
                      signif.SG=rep(NA, n.temp), signif.MR=rep(NA, n.temp))
simres$N.A <- kronecker(N.A, rep(1, n.B))
simres$N.B <- rep(N.B, n.A)

#### Run the scenario
  for(i in 1:n.temp)
  {
  #### Run the scenario
  N.A.temp <- simres$N.A[i]
  N.B.temp <- simres$N.B[i]
  temp <- sim.meta(N.sim, N.A.temp, N.B.temp, RR.A, RR.B, SE.A, SE.B, T.A, T.B, V.A, V.B)
    
  #### Save the results
  simres[i, 3:16] <- as.numeric(t(temp))[c(1:12, 14:15)]
  }


#### View and save the results
simres
write.csv(simres, file="Simulation results.csv")
