### Script: Bayesian regression linear accounting for measurement error
### Plus uncertainity in site level (spatial random effect = eps.site)
### Based on Kery & Royle (Applied hierarchical modeling in ecology: Analysis of distribution, abundance and species richness in R and BUGS Volume 1))
### and Maxwell Joseph (https://mbjoseph.github.io/posts/2018-12-23-bayesian-model-ii-regression/)
### Code by Ardiantiono and Nicolas J. Deere
### -----------------------------------------------------------------------

# Clear workspace
rm(list=ls())

# Load required library
library("jagsUI")

# Set working directory (where the all data required to run the models are stored)
setwd("~") #Set up file path
getwd()

# Load data
data <- read.csv("Species_Biodiversity.csv", header=TRUE)
names(data)

##########################################################################
## Regression linear accounting for response & covariate measurement error
## For Community Occupancy and Species Richness calculation
## As they have standard deviation from occupancy model output 
## -----------------------------------------------------------------

# Extract estimates of N or biodiversity parameter (this is the response variable)
N.pm1 <- data[c(9, 11)] #Posterior mean of biodiversity parameters
N.psd1 <- data[c(10, 12)] #Posterior SD of biodiversity parameters
N.list1 <- c("Community occupancy", "Species richness")

# Extract estimates of SP or species occupancy  (this is the covariate variable)
SP.pm <- data[,c(3,5,7)] # Posterior means of umbrella species occupancy
SP.psd <- data[,c(4,6,8)] # Posterior SD's of umbrella species occupancy
species.list <- c("SP1", "SP2", "SP3")

# Create a sequence for prediction
pred.occu <- seq(0, 1, 0.01)  #no need to standardize


#Create template for df 
beta_cat1 <- data.frame(matrix(ncol = 14, nrow = 0))

####################################
for(j in 1:length(N.list1)){
  for(k in 1:length(species.list)) {
    
    str(win.data <- list(N = N.pm1[,j], psd = N.psd1[,j], n = length(N.pm1[,j]),
                         SP = SP.pm[,k], ssd = SP.psd[,k], sp.err = mean(SP.psd[,k]),           # sp.err = mean standard deviation of sp occupancy covariate. 
                         pred.occu = pred.occu, npred = length(pred.occu)))
    
# Define model in BUGS language 
sink("meta.analysis.txt") 
cat("
model{

# Priors 
for(v in 1:3){ # Priors for intercept and polynomial coefficients 
  beta[v] ~ dnorm(0, 0.0001)
  }
tau.site <- pow(sd.site, -2) 
sd.site ~ dunif(0,100)


# Likelihood 
for (i in 1:n){ 
    truex[i] ~ dnorm(0, sp.err)
    SP[i] ~ dnorm(truex[i], tau.ssd[i])
    tau.ssd[i] <- pow(ssd[i], -2)
    
    N[i] ~ dnorm(muN[i], tau.psd[i])
    tau.psd[i] <- pow(psd[i], -2)  
    
    #-----------------------------------------------------------------------------------------------
    muN[i] <- beta[1] + beta[2]*truex[i] + beta[3]*pow(truex[i],2)  + eps.site[i] #Remove eps.site (accounting for uncertainity in site)
    eps.site[i] ~ dnorm(0, tau.site)
    
    #--------------------------------------------------------------------------------
        # Assess model fit using a sums-of-squares-type discrepancy
        residual[i] <- N[i]-muN[i]             # Residuals for observed data
        predicted[i] <- muN[i]                # Predicted values
        sq[i] <- pow(residual[i], 2)     # Squared residuals for observed data
    
        # Generate replicate data and compute fit stats for them
        y.new[i] ~ dnorm(muN[i], tau.psd[i]) # one new data set at each MCMC iteration
        sq.new[i] <- pow(y.new[i]-predicted[i], 2)   # Squared residuals for new data
}
  
 fit <- sum(sq[])                       # Sum of squared residuals for actual data set
 fit.new <- sum(sq.new[])             # Sum of squared residuals for new data set
 test <- step(fit.new - fit)               # Test whether new data set more extreme
 bpvalue <- mean(test)                    # Bayesian p-value

# Get predictions for plot 
for(i in 1:npred){ 
  Npred[i] <- beta[1] + beta[2]*pred.occu[i] +beta[3]*pow(pred.occu[i],2) 
  } 

} # model end-------------------------------------------------------------------------
",fill=TRUE) 
    sink()    
    
    # Initial values, params monitored, and MCMC settings 
    inits <- function() list(beta = rnorm(3), sd.site = runif(1)) 
    params <- c("beta", "sd.site", "Npred", "fit", "fit.new", "bpvalue") 
    #ni <- 100 ; nt <- 10 ; nb <- 10 ; nc <- 3
    ni <- 50000 ; nt <- 25 ; nb <- 25000 ; nc <- 3
    
    # Call JAGS and summarize posterior 
    out <- jags(win.data, inits, params, "meta.analysis.txt", 
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb) 
    
    jags.summ <- out$summary
    
    #Combine the output
    biodiv <- N.list1[j]
    sp.name <- species.list[k]
    
    beta.est1 <- as.data.frame(cbind(jags.summ, biodiv, sp.name))
    beta_cat1 <- rbind(beta_cat1, beta.est1)
  }
}

write.csv(beta_cat1, file=paste("XXX", ".csv", sep="")) #add file name
save(out, file=paste("XXX", ".RData", sep="")) #add file name



####################################################################
####################################################################
## Regression linear accounting for only covariate measurement error
## For FD and PD calculation
## -----------------------------------------------------------------

# Extract estimates of N or biodiversity parameter (this is the response variable)
N.pm <- data[c(13,14)] #Estimate of biodiversity parameters
N.list <- c("SES FDis", "SES MPD")

# Extract estimates of SP or species occupancy  (this is the covariate variable)
SP.pm <- data[,c(3,5,7)] # Posterior means of umbrella species occupancy
SP.psd <- data[,c(4,6,8)] # Posterior SD's of umbrella species occupancy
species.list <- c("SP1", "SP2", "SP3")

# Create a sequence for prediction
pred.occu <- seq(0, 1, 0.01)  #no need to standardize

#Create template for df 
beta_cat <- data.frame(matrix(ncol = 14, nrow = 0))

####################################
for(j in 1:length(N.list)){
  for(k in 1:length(species.list)) {
    
    str(win.data <- list(N = N.pm[,j], n = length(N.pm[,j]),
                         SP = SP.pm[,k], ssd = SP.psd[,k], sp.err = mean(SP.psd[,k]),           # sp.err = mean standard deviation of sp occupancy covariate. 
                         pred.occu = pred.occu, npred = length(pred.occu)))
    
# Define model in BUGS language 
sink("meta.analysis.txt") 
cat("
model{

# Priors 
for(v in 1:3){ # Priors for intercept and polynomial coefficients 
  beta[v] ~ dnorm(0, 0.0001)
  }
tau.site <- pow(sd.site, -2) 
sd.site ~ dunif(0,100)

# In Maxwell they set up prior for tau.ssd ~ dunif and tau.psd <- pow here... 
# So only one value for tau, not in each site

# Likelihood 
for (i in 1:n){ 
    truex[i] ~ dnorm(0, sp.err)
    SP[i] ~ dnorm(truex[i], tau.ssd[i])
    tau.ssd[i] <- pow(ssd[i], -2)
    
    N[i] ~ dnorm(muN[i], tau.site) 
    
    #-----------------------------------------------------------------------------------------------
    muN[i] <- beta[1] + beta[2]*truex[i] + beta[3]*pow(truex[i],2)  + eps.site[i] #Remove eps.site (accounting for uncertainity in site because already using site random effects)
    eps.site[i] ~ dnorm(0, tau.site)
    
    #--------------------------------------------------------------------------------
        # Assess model fit using a sums-of-squares-type discrepancy
        residual[i] <- N[i]-muN[i]             # Residuals for observed data
        predicted[i] <- muN[i]                # Predicted values
        sq[i] <- pow(residual[i], 2)     # Squared residuals for observed data
    
        # Generate replicate data and compute fit stats for them
        y.new[i] ~ dnorm(muN[i], tau.site) # one new data set at each MCMC iteration
        sq.new[i] <- pow(y.new[i]-predicted[i], 2)   # Squared residuals for new data
}
  
 fit <- sum(sq[])                       # Sum of squared residuals for actual data set
 fit.new <- sum(sq.new[])             # Sum of squared residuals for new data set
 test <- step(fit.new - fit)               # Test whether new data set more extreme
 bpvalue <- mean(test)                    # Bayesian p-value

# Get predictions for plot 
for(i in 1:npred){ 
  Npred[i] <- beta[1] + beta[2]*pred.occu[i] +beta[3]*pow(pred.occu[i],2) 
  } 

} # model end-------------------------------------------------------------------------
",fill=TRUE) 
    sink()    
    
    
    # Initial values, params monitored, and MCMC settings 
    inits <- function() list(beta = rnorm(3), sd.site = runif(1)) 
    params <- c("beta", "sd.site", "Npred", "fit", "fit.new", "bpvalue") 
    #ni <- 100 ; nt <- 10 ; nb <- 10 ; nc <- 3
    ni <- 50000 ; nt <- 25 ; nb <- 25000 ; nc <- 3
    
    # Call JAGS and summarize posterior 
    out <- jags(win.data, inits, params, "meta.analysis.txt", 
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb) 
    
    jags.summ <- out$summary
    
    #Combine the output
    biodiv <- N.list[j]
    sp.name <- species.list[k]
    
    beta.est <- as.data.frame(cbind(jags.summ, biodiv, sp.name))
    beta_cat <- rbind(beta_cat, beta.est)
  }
}

write.csv(beta_cat, file=paste("XXX", ".csv", sep="")) #add file name
save(out, file=paste("XXX", ".RData", sep="")) #add file name



####################
#Plot the prediction
plot(SP.pm, N.pm, xlab = "Species occupancy", 
     ylab = "Estimated community occupancy", ylim = c(0, max(N.pm)), xlim=c(0,max(SP.pm)), frame = F) 
lines(seq(0, 1, 0.01), out$mean$Npred, col = "blue", lwd = 3)                             # Line predicted by model
matlines(seq(0, 1, 0.01), out$summary[5:105,c(3, 7)], col = "blue", lwd = 2, 
         lty= "dashed")
lines(smooth.spline(N.pm ~ SP.pm), col = "red", lwd = 2)                   # Standard regression line

help(lines)



