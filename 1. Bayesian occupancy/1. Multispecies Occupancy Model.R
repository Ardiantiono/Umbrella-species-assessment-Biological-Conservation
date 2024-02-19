### R code detailing data processing for hierarchical occupancy model to explore ecological 
### correlates of medium/large terrestrial mammal occurrence 
### With Data Augmentation & Spatial random effect
### BUGS specification of the hierarchical model framework is included, to be implemented in
### JAGS, operated via the R console using the "jagsUI" library
### Code produced by Nicolas J. Deere adapted by Ardiantiono
#==================================================================================================================
# 1a. Data Preparation:
#==================================================================================================================
rm(list=ls())

# Set working directory (where the all data required to run the models are stored)
setwd("~") #fill with directory path
getwd()

# Read and check occupancy detection data
# Note that this dataframe has had all NA occasions replaced with zeroes, these must be
# added manually. 
# Incorporating NAs in the model at this point will cause R to treat detections as categorical
# with three levels ("0", "1", "NA").
det.mat <- read.table("Combined_Detection_AllSpecies_NoNAs.csv", header=TRUE, sep=",", na.strings=TRUE)
names(det.mat)
head(det.mat)
str(det.mat)

# Read in dataframe containing detection data including NAs - we will use these data to
# extract and add NA values at specific sites
# Note that we subset the data for a single species as this will be multiplied across the
# third dimension of the array (relating to species)
df.nas <- read.table("Combined_Detection_AllSpecies.csv", header=TRUE, sep=",", na.strings=TRUE,
                     colClasses = c("factor", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", "numeric", "numeric",  "factor")) #specify class perhaps cause NAs
df.nas <- df.nas[1:48,1:5] #use only SP1 data row
df.nas <- as.matrix(df.nas[,2:19])
head(df.nas)
str(df.nas)


#Calculate number of sightings for each species
total.count <- tapply(det.mat$Total.Detection, det.mat$Species, sum) #start from number then character
total.count

#Find the number of unique species
uspecies <- as.character(unique(det.mat$Species)) #retain the character format
uspecies

#n is the number of observed species
n <- length(uspecies)

#Augment dataset: create number of undetected species
nz <- 23 #undetected number (species > 1 kg from Sibarani et al (2019))
M <- n+nz #superpopulation number (28 sp)

#Find the number of unique sampling locations
upoints = as.character(unique(det.mat$Site_ID))
#J is the number of sampled points
J=length(upoints)

# Convert dataframe into a 3-dimensional array
# Site (row) * Temporal Replicate (column) * Species (3rd dimension)
# Dimensions renamed for consistency and ease of recognition

y <- simplify2array(by(det.mat[2:19], det.mat$Species, as.matrix, row.names=det.mat$Site_ID))
dimnames(y)[[1]] <- upoints #name dim 1
occasions <- rep(1:18,1) #create vector
nrep <- length(occasions)
dimnames(y)[[2]] <- occasions #name dim 2
y #see the result

#Create matrix with superpopulation
yaug <- array(0, dim=c(J, nrep, n+nz)) #18 = number of SO
yaug[,,1:n] <- y

# Add in the missing replicates with NAs derived from the df.nas matrix
for (i in 1: dim(yaug)[3]) {               #[3] means the detection matrix, not the header
  yaug[,,i][is.na(df.nas)] <- NA
}
summary(yaug)

# Now we create a matrix of all ones adding the NAs from failed sampling occasions
# This will be used to sum the number of sampling occasions per site (K)
y.one = matrix(1, nrow=48, ncol=18)
y.one[is.na(df.nas)] <- NA
K=apply(y.one,1,sum, na.rm=TRUE) # creates a vector, "K", containing number of detections (1s) per row
K=as.vector(K) # converts to class vector

# 1b. Covariate processing
#==================================================================================================================
# Read and check covariates dataframe
covars <- read.table("OM_Covariates.csv", header=TRUE,sep=",",na.strings=c("NA"))
names(covars)
str(covars)

# Scale continuous predictor variables
# Write custom centre and scaling function
var.stdz <- function(x) {
  m.x <- mean(x)
  sd.x <- sd(x)          #same with scale function?
  (x-m.x)/sd.x
}

covars$S_trap_effort <- var.stdz(covars$Trap_Effort)

covars$S_tri_1500 <- var.stdz(covars$tri_1500)
covars$S_elev_1000 <- var.stdz(covars$elev_1000)

covars$S_biomass_1000 <- var.stdz(covars$biomass_1000)
covars$S_prop_forest_1500 <- var.stdz(covars$prop_forest_1500)
covars$S_d_water <- var.stdz(covars$d_water)

covars$S_ttcsm <- var.stdz(covars$ttcsm)

# Remove non-standardized variables from dataframe
names(covars)
covars <- covars[,-c(7,9:14)] #double clean non std columns
names(covars)

# Create variables for unique number of sampling blocks contained within the sampling design
n.block <- length(unique(covars$fSite))

# Create vectors linking each site to a specific year category (no need)
block.counter <- covars$fSite

# Create variables for unique number of CT type used 
# 1: Panthera IV, 2: Panthera V, 3: Reconyx, 4: Bushnell
l.type <- length(unique(covars$CT_type))

# Create a list and grouped covariates (only occupancy covariates)
names(covars)
covar.list <-  covars[,9:14] 

# Create new covariate dataframe to use for the graphing procedure
vars <- read.table("OM_Covariates_FINAL.csv", header=TRUE,sep=",",na.strings=c("NA"))
names(vars)
var.list <- vars[,9:14] 


# 1c. Bayesian Hierarchical Modelling code (written in BUGS language)
# Edited for global occu model
#==================================================================================================================
# Set working directory to analysis folder
# Dorazio-Royle community occupancy model (with data augmentation)
# Half-cauchy priors implemented as weakly informative priors to shrink parameter estimates 
# towards zero and mitigate variance inflation (following Broms et al. 2016; Ecology) 
# Add = CT type as effect of detectability (use nested indexing; type as intercept factor)

#some notes:
#mu (mean), tau (variance/precision)
#dunif() = uniform distribution; min to max range
#pow () = power transformation for precision
#dbern = #bernoulli distribution
#step = stepwise
#sigma = variance

sink("OM_AllSpecies_multivariate.txt")
cat("
    model{
    
    #Priors
    omega ~ dunif(0,1)
    
     # Priors for species-specific effects on occupancy and detection
    #===============================================================
    #for alpha (occupancy)
    for(i in 1:(n.sp+n.zeroes)){
       for(s in 1:n.PsiParams){
           alpha[i,s] ~ dnorm(mu.alpha[s], tau.alpha[s])
       }
    }
    
    #for beta intercept
    for(i in 1:(n.sp+n.zeroes)){
        for(s in 1:n.type){
              beta.int[s,i] ~ dnorm(mu.beta.int[s], tau.beta.int[s])
          }
       }

    #for beta 
    for(i in 1:(n.sp+n.zeroes)){
        for(s in 1:n.DetParams){
            beta[i,s] ~ dnorm(mu.beta[s], tau.beta[s])
       }
    }
    
        
    # Hyperpriors for process/occupancy model
    #========================================
    for(s in 1:n.PsiParams){
        mu.alpha[s] ~ dnorm(0, 0.01) #check this if rhat high (uniform?)
        sigma.alpha[s] ~ dunif(0, 10)
        tau.alpha[s] <- pow(sigma.alpha[s], -2)
    }

    # Hyperpriors for observation/detection model
    #============================================
    for(s in 1:n.DetParams){
        mu.beta[s] ~ dnorm(0, 0.01)
        sigma.beta[s] ~ dunif(0, 10)
        tau.beta[s] <- pow(sigma.beta[s], -2)
    }
    
    #beta intercept
      for(s in 1:n.type){
        mu.beta.int[s] ~ dnorm(0, 0.01)
        sigma.beta.int[s] ~ dunif(0, 10)
        tau.beta.int[s] <- pow(sigma.beta.int[s], -2)
    }


    # Hyperpriors/priors for spatial random effects
    #===========================================================
    for(i in 1:(n.sp+n.zeroes)) { 
    # Random block effects for occupancy
        for(block in 1:n.block){
            delta[block, i] ~ dnorm(0, delta.tau[i])
        }
      delta.tau[i] <- pow(delta.sd[i], -2)
      delta.sd[i] ~ dunif(0, 10)
    }
   
    # Superpopulation process: Ntotal species sampled out of M available
    #===========================================================
    for(i in 1:(n.sp+n.zeroes)) {
    w[i] ~ dbern(omega)
    }
    
    # Ecological model for true occurence (process model; I used 6 covariates; tri, ttcsm, and elev use quadratic)
    #====================================================
    for(i in 1:(n.sp+n.zeroes)) {
        for(j in 1:n.sites) {
            logit(psi[j,i]) <- alpha[i,1] + alpha[i,2]*cov1.psi[j] + alpha[i,3]*cov1.psi2[j] +
                               alpha[i,4]*cov2.psi[j] + alpha[i,5]*cov2.psi2[j] +
                               alpha[i,6]*cov3.psi[j] + alpha[i,7]*cov4.psi[j] +
                               alpha[i,8]*cov5.psi[j] + alpha[i,9]*cov6.psi[j] + alpha[i,10]*cov6.psi2[j] +
                               delta[block.counter[j], i]
            mu.psi[j,i] <- w[i]*psi[j,i] #update: recalculate mu.psi with w 
            z[j,i] ~ dbern(mu.psi[j,i])
        }
    }
    
    # Observation model for replicated detection/non-detection observations (trap effort + CT type)
    #======================================================================
    for(i in 1:(n.sp+n.zeroes)) {
        for(j in 1:n.sites) {    
            for(k in 1:n.reps[j]) {
            logit(p[j,k,i]) <- beta.int[cov0.p[j],i] + beta[i,1]*cov1.p[j] + beta[i,2]*cov2.p[j] #issue here, beta.int unknown!
            mu.p[j,k,i] <- p[j,k,i]*z[j,i]
            y[j,k,i] ~ dbern(mu.p[j,k,i]) 
            }
        }
    }
    
    # Calculate Pearson's chi-squared residuals to assess goodness of fit
    # Based on Kery and Royle: Applied hierarchical modelling in ecology, pp. 235
    # Calculate the observed and expected residuals
    # Add small value to prevent division by zero
    #============================================
    for(i in 1:(n.sp+n.zeroes)) {
        for(j in 1:n.sites) { 
            for(k in 1:n.reps[j]) {
        
            #calculate for each SO (split the formula into 2 steps)
            y.sim[j,k,i] ~ dbern(mu.p[j,k,i])                                                # Replicated dataset
            chi.obs[j,k,i] <- abs(y[j,k,i] - mu.p[j,k,i])        # Observed dataset. Add small value to denominator to prevent division by zero
            chi.sim[j,k,i] <- abs(y.sim[j,k,i] - mu.p[j,k,i])    # Expected dataset. Add small value to denominator to prevent division by zero
            chi2.obs[j,k,i] <- pow(chi.obs[j,k,i], 2)
            chi2.sim[j,k,i] <- pow(chi.sim[j,k,i], 2)
    
            }
            
            
         #calculate for each site
        chi2.obs.sum[j,i] <- sum(chi2.obs[j,1:n.reps[j],i])
        chi2.sim.sum[j,i] <- sum(chi2.sim[j,1:n.reps[j],i])  
          }
            
    # Calculate chi-squared discrepancy for each species 
    #===================================================
      fit.sp.obs[i] <- sum(chi2.obs.sum[,i])                                     # Species-specific fit statistic for actual dataset
      fit.sp.sim[i] <- sum(chi2.sim.sum[,i])                                     # Species-specific fit statistic for simulated dataset
      c.hat.sp[i] <- fit.sp.obs[i]/fit.sp.sim[i]
      bpv.sp[i] <- step(fit.sp.sim[i] - fit.sp.obs[i])
    }
    
    
    # Calculate overall chi-squared discrepency measure
    # Alt. to bpv after model output run mean(out$sims.list$fit.obs > out$sims.list$fit.sim)
    #=======================================================================================
    fit.obs <- sum(chi2.obs.sum[1:n.sites, 1:n.sp]) #using n.species only
    fit.sim <- sum(chi2.sim.sum[1:n.sites, 1:n.sp])
    c.hat <- fit.obs/fit.sim
    bpv <- step(fit.sim - fit.obs)
    
    # Derived quantities
    # Number of occupied sites
    #=========================
    for(i in 1:(n.sp+n.zeroes)) {
        Nocc.fs[i] <- sum(z[,i])
        }
    
    # Number of species occurring at each site
    #=========================================
    for(j in 1:n.sites) {
        Nsite[j] <- inprod(z[j,1:(n.sp+n.zeroes)],w[1:(n.sp+n.zeroes)])
    }
    
    # Number of unseen species and metacommunity size size
    #=========================================
    n0 <- sum(w[(n.sp+1):(n.sp+n.zeroes)])
    N <- n.sp + n0
        
    }
    ",fill=TRUE)
sink()


# 1d. Model preparation and initialization
#==================================================================================================================
# Load library to call JAGS from R
library(jagsUI)
library(plyr)
library(reshape2)
library(boot)  # for logit functions
library(verification)  # to calc AUC statistic (used in functions)

# Specify number of covariates for occupancy and detection models
n.PsiParams <- 10 #intercept, 3 covariates with quadratic (3 x 2), 3 covariates with linear (3 x 1)
n.DetParams <- 2 #not include intercept (CT type); trap_effort, and elevation
n.type <- 4


# Specify the fundamentals of the MCMC chains
#ni <- 100000; nb <- 50000; nt <- 50; nc <- 3 # Final iterations
ni <- 1000; nb <- 100; nt <- 10; nc <- 3 # Exploratory iterations (nt=thinning rate, nc= chain)

# Create a reformatted version of the occupancy matrix (n.sp x n.sites x n.reps)
Y <- aperm(yaug, c(3,1,2))

# Generate probability expressions (if we call WAIC script to calculate wAIC)
exp.tmp <- createExpressions(mod.type="AllData")
expressPsi <- exp.tmp$expressPsi
expressDet <- exp.tmp$expressDet

# For reproduceability
set.seed(2000)

# Run univariate models across all covariates 
{
  m.fit <- setNames(data.frame(matrix(ncol = 8, nrow = 1)), c("Model", "DIC", "BPV", "Chat", "lppd", "pD", "WAIC", "CPO"))
  ### Analyze the data
  ###=================
  # Specify the parameters to be monitored  
  params <- c("omega", "mu.alpha", "mu.beta", "mu.beta.int", "sigma.alpha", "sigma.beta", "sigma.beta.int",
              "alpha", "beta.int", "beta", "z", "delta", "delta.tau", "Nocc.fs", 
              "Nsite",  "w", "n0", "N", "fit.obs", "fit.sim", "c.hat", "bpv.sp", "bpv")
  
  # Create standardized numeric character strings for the covariates that change within the loop
  cov1.psi <- as.vector(covars$S_tri_1500)
  cov1.psi2 <- cov1.psi*cov1.psi
  
  cov2.psi <- as.vector(covars$S_elev_1000)
  cov2.psi2 <- cov2.psi*cov2.psi
  
  cov3.psi <- as.vector(covars$S_biomass_1000)
  
  cov4.psi <- as.vector(covars$S_prop_forest_1500)
  
  cov5.psi <- as.vector(covars$S_d_water)
  
  cov6.psi <- as.vector(covars$S_ttcsm)
  cov6.psi2 <- cov6.psi*cov6.psi
  
  cov0.p <- as.vector(covars$CT_type)
  cov1.p <- as.vector(covars$S_trap_effort)
  cov2.p <- as.vector(covars$S_elev_1000) #no need quadratic term
  
  
  # Load the data
  data <- list(n.sp=n, n.zeroes=nz, n.sites=J, n.reps=K, y=yaug, n.PsiParams=n.PsiParams, 
               n.type=n.type, n.DetParams=n.DetParams,  cov1.psi= cov1.psi, 
               cov1.psi2 = cov1.psi2, cov2.psi= cov2.psi, cov2.psi2= cov2.psi2,
               cov3.psi= cov3.psi, cov4.psi= cov4.psi, cov5.psi= cov5.psi, cov6.psi= cov6.psi,
               cov6.psi2= cov6.psi2, cov0.p=cov0.p, cov1.p=cov1.p, cov2.p=cov2.p, n.block=n.block,
               block.counter=block.counter)  
  # Specify initial values
  wst <- rep(1, M)
  zst <- array(1, dim=c(J, M))
  inits <- function() list(z=zst, w=wst, mu.alpha=rnorm(n.PsiParams), sigma.alpha=runif(n.PsiParams, 0, 10),
                           mu.beta=rnorm(n.DetParams), mu.beta.int=rnorm(n.type), sigma.beta=runif(n.DetParams, 0, 10),
                           sigma.beta.int=runif(n.type, 0, 10), delta.sd=runif(M,0,10)) 
  #Run the model and extract the results
  print(paste('HMSOM Univariate Modelling: Mammals', 'Time: ', Sys.time()))
  out <- jags(data, inits, params, "OM_AllSpecies_multivariate.txt", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE)
  save(out, file=paste("HMSOM_AllSpecies_DA", ".RData", sep=""))
  jags.summ <- out$summary
  write.csv(jags.summ, file=paste("HMSOM_AllSpecies_DA", ".csv", sep=""))
}
  