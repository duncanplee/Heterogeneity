sim.meta <- function(N.sim, N.A, N.B, RR.A, RR.B, SE.A, SE.B, t.A, t.B, nu.A, nu.B)
{  
#### Compute the log relative risks
beta.A <- log(RR.A)
beta.B <- log(RR.B)


#### Create matrices to store the results
results <- array(NA, c(N.sim, 11))
colnames(results) <- c("est.global", "est.SG", "est.MR", 
                       "LCI.global", "UCI.global", "LCI.SG", "UCI.SG", "LCI.MR", "UCI.MR",
                       "p-value.SG", "p-value.MR")
results <- as.data.frame(results)



#############################
#### Run the simulation study
#############################
  for(i in 1:N.sim)
  {
  #### Generate data for a meta analysis
  estimates.A <- rnorm(n=N.A, mean=beta.A, sd=t.A)  
  estimates.se.A <- rnorm(n=N.A, mean=SE.A, sd=nu.A)
  estimates.se.A[estimates.se.A<0] <- 0.0001
  
  estimates.B <- rnorm(n=N.B, mean=beta.B, sd=t.B)  
  estimates.se.B <- rnorm(n=N.B, mean=SE.B, sd=nu.B)
  estimates.se.B[estimates.se.B<0] <- 0.0001

  dat <- data.frame(Region=c(rep("A", N.A), rep("B", N.B)), 
                    estimate=c(estimates.A, estimates.B),
                    se=c(estimates.se.A, estimates.se.B))  

  #### Global random effects model to all studies
  mod1 <- metagen(TE=estimate, seTE=se, data=dat, sm="RR", comb.random=TRUE, method.tau="DL", hakn=TRUE, backtransf=TRUE, prediction=TRUE, level.predict=0.95)
  results[i, 1] <- exp(mod1$TE.random)
  results[i, 4:5] <- exp(c(mod1$lower.random, mod1$upper.random))


  #### Sub-group analysis by region
  mod2 <- update.meta(mod1, data=dat, subgroup = Region, tau.common = FALSE)
  results[i, 2] <- exp(mod2$TE.random.w[1])
  results[i, 10] <- mod2$pval.Q.b.random
  results[i, 6:7] <- c(exp(mod2$lower.random.w[1]), exp(mod2$upper.random.w[1]))

  
  #### Meta-regression with region as a covariate - effect estimates and CIs
  mod3a <- metareg(x=mod1, formula=~factor(Region), intercept=FALSE, hakn=TRUE)
  results[i, 3] <- exp(mod3a$beta[1])
  results[i, 8:9] <- exp(c(mod3a$ci.lb[1], mod3a$ci.ub[1]))

  
  #### Meta-regression for significance test
  mod3b <- metareg(x=mod1, formula=~factor(Region), intercept=TRUE, hakn=TRUE)
  results[i, 11] <- mod3b$QMp
  }


#####################################
#### Summarise and return the results
#####################################
#### Set up a results matrix
results.final <- array(NA, dim=c(5,3)) 
colnames(results.final) <- c("Global", "SG", "MR")
rownames(results.final) <- c("Bias", "RMSE", "Coverage", "Width", "Signif")


#### Save the results
results.final[1, ] <- 100 * apply(log(results[ ,1:3]) - beta.A , 2, mean) / beta.A
results.final[2, ] <- 100 * sqrt(apply((log(results[ ,1:3]) - beta.A)^2 , 2, mean)) / beta.A
results.final[3, 1] <- 100 * sum(results$LCI.global < RR.A & results$UCI.global > RR.A) / N.sim
results.final[3, 2] <- 100 * sum(results$LCI.SG < RR.A & results$UCI.SG > RR.A) / N.sim
results.final[3, 3] <- 100 * sum(results$LCI.MR < RR.A & results$UCI.MR > RR.A) / N.sim
results.final[4, 1] <- mean(results$UCI.global - results$LCI.global)
results.final[4, 2] <- mean(results$UCI.SG - results$LCI.SG)
results.final[4, 3] <- mean(results$UCI.MR - results$LCI.MR)
results.final[5, 2] <- mean(results$`p-value.SG` < 0.05)
results.final[5, 3] <- mean(results$`p-value.MR` < 0.05)
return(results.final)
}

