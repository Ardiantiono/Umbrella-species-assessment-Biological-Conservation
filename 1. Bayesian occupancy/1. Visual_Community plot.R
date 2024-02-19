### Functions to generate occupancy graphs
### Code developed by Nicolas J. Deere adapted by Ardiantiono

rm(list=ls())

# Set working directory (where the all data required to run the models are stored)
setwd("~") #Set up file path
getwd()

jags.out <- load("~//HMSOM_AllSpecies_DA.RData") #Call the model output 
names(out) #it reads as out not jags.out...

######################################

library(plyr) #for function round.any

# call covariate ##########
vars <- read.table("OM_Covariates.csv", header=TRUE,sep=",",na.strings=c("NA"))
names(vars)
var.list <- vars[,8:13] 

psi.covar1 <- as.vector(var.list$tri_1500)
psi.covar2 <- psi.covar1*psi.covar1
psi.covar3 <- as.vector(var.list$elev_1000)
psi.covar4 <- psi.covar3*psi.covar3
psi.covar5 <- as.vector(var.list$biomass_1000)
psi.covar6 <- as.vector(var.list$prop_forest_1500)
psi.covar7 <- as.vector(var.list$d_water)
psi.covar8 <- as.vector(var.list$ttcsm)
psi.covar9 <- psi.covar8*psi.covar8


# TRI #######
{
  # Create objects to define the parameters of the model prediction
  covar.min.value <- round_any(min(psi.covar1), 1, floor)
  covar.max.value <- round_any(max(psi.covar1), 1, ceiling)
  
  # Calculate covariate values for prediction
  var.pred <- seq(covar.min.value, covar.max.value,,500)
  var.mean <- mean(psi.covar1)
  var.sd <- sd(psi.covar1)
  var.pred.stdz <- (var.pred - var.mean)/var.sd
  
  # Create a temp folder storing species-specific posterior distributions of all monitored parameters
  # Compute predictions based on the number of MCMC samples
  n.samp <- out$mcmc.info$n.samples
  pred <- array(NA, dim=c(500, n.samp,1))
  for(kk in 1:n.samp){
    pred[,kk,1] <- plogis(out$sims.list$mu.alpha[,1][kk] + out$sims.list$mu.alpha[,2][kk] * var.pred.stdz
                   + out$sims.list$mu.alpha[,3][kk] * var.pred.stdz^2)
    }
  
  #also create vector of occupancy estimate in each station
  site.pred.stdz <- (psi.covar1-var.mean)/var.sd
  
  occu_pred <- plogis(out$mean$mu.alpha[1] + out$mean$mu.alpha[2] * site.pred.stdz 
                      + out$mean$mu.alpha[3] * site.pred.stdz^2)
  y.jitter <- jitter(occu_pred, factor=1000) #add noise, so points not overlap
  #the value really depends on the variance of data, so need to calibrate
  
  #then join the  prediction with occu_pred
  site_occu <- cbind(psi.covar1, occu_pred)
  names(site_occu)
  colnames(site_occu) <- c("Covariate", "Psi")
  site_occudf <- as.data.frame(site_occu)
  
  
  #create the plot!
  pmC <- apply(pred, c(1,3), mean)
  criC <- apply(pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
  
  # Plot posterior mean and a random sample of 100 from posterior of regression
  selection <- sample(1:n.samp,500)
  #selection <- sample(1:n.samp,10) #exploratory
  jpeg(file=paste("~", #Directory path
                  "Community Plot TRI", ".jpeg", sep=""), width=1000, height=800)
  par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,1.1))
  matplot(var.pred, pred[,selection,], ylab="Mean Occupancy", xlab="TRI", 
          type = "l", lty = 1, lwd=1, col="grey", ylim=c(0,1), xlim=c(covar.min.value, covar.max.value),
          cex.axis=1.5, cex.lab=2)
  lines(var.pred, apply(pred[,,],1,mean), lwd=3, col="blue")
  matlines(var.pred, t(criC[,,]), col="blue", lty=3, lwd=2)
  points(x = site_occudf$Covariate, y = y.jitter, pch = 16,
         col = "black")
  mtext(paste("Community-level Mammalian Occupancy: ", "TRI", sep=" "), side=3, line=-2.5, cex=1.5, outer=TRUE)
  dev.off()
}


### Elevation ###############
{
  # Create objects to define the parameters of the model prediction
  covar.min.value <- round_any(min(psi.covar3), 1, floor)
  covar.max.value <- round_any(max(psi.covar3), 1, ceiling)
  
  # Calculate covariate values for prediction
  var.pred <- seq(covar.min.value, covar.max.value,,500)
  var.mean <- mean(psi.covar3)
  var.sd <- sd(psi.covar3)
  var.pred.stdz <- (var.pred - var.mean)/var.sd
  
  # Create a temp folder storing species-specific posterior distributions of all monitored parameters
  # Compute predictions based on the number of MCMC samples
  n.samp <- out$mcmc.info$n.samples
  pred <- array(NA, dim=c(500, n.samp,1))
  for(kk in 1:n.samp){
    pred[,kk,1] <- plogis(out$sims.list$mu.alpha[,1][kk] + out$sims.list$mu.alpha[,4][kk] * var.pred.stdz
                          + out$sims.list$mu.alpha[,5][kk] * var.pred.stdz^2)
  }
  
  #also create vector of occupancy estimate in each station
  site.pred.stdz <- (psi.covar3-var.mean)/var.sd
  
  occu_pred <- plogis(out$mean$mu.alpha[1] + out$mean$mu.alpha[4] * site.pred.stdz 
                      + out$mean$mu.alpha[5] * site.pred.stdz^2)
  
  y.jitter <- jitter(occu_pred, factor=1000) #add noise, so points not overlap
  #the value really depends on the variance of data, so need to calibrate
  
  #then join the  prediction with occu_pred
  site_occu <- cbind(psi.covar3, occu_pred)
  names(site_occu)
  colnames(site_occu) <- c("Covariate", "Psi")
  site_occudf <- as.data.frame(site_occu)
  
  
  #create the plot!
  pmC <- apply(pred, c(1,3), mean)
  criC <- apply(pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
  
  # Plot posterior mean and a random sample of 100 from posterior of regression
  selection <- sample(1:n.samp,500)
  #selection <- sample(1:n.samp,10) #exploratory
  jpeg(file=paste("~", #Directory path
                  "Community Plot Elev", ".jpeg", sep=""), width=1000, height=800)
  par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,1.1))
  matplot(var.pred, pred[,selection,], ylab="Mean Occupancy", xlab="Elevation", 
          type = "l", lty = 1, lwd=1, col="grey", ylim=c(0,1), xlim=c(covar.min.value, covar.max.value),
          cex.axis=1.5, cex.lab=2)
  lines(var.pred, apply(pred[,,],1,mean), lwd=3, col="blue")
  matlines(var.pred, t(criC[,,]), col="blue", lty=3, lwd=2)
  points(x = site_occudf$Covariate, y = y.jitter, pch = 16,
         col = "black")
  mtext(paste("Community-level Mammalian Occupancy: ", "Elevation", sep=" "), side=3, line=-2.5, cex=1.5, outer=TRUE)
  dev.off()
}

### Biomass ####
 {
  # Create objects to define the parameters of the model prediction
  covar.min.value <- round_any(min(psi.covar5), 1, floor)
  covar.max.value <- round_any(max(psi.covar5), 1, ceiling)
  
  # Calculate covariate values for prediction
  var.pred <- seq(covar.min.value, covar.max.value,,500)
  var.mean <- mean(psi.covar5)
  var.sd <- sd(psi.covar5)
  var.pred.stdz <- (var.pred - var.mean)/var.sd
  
  # Create a temp folder storing species-specific posterior distributions of all monitored parameters
  # Compute predictions based on the number of MCMC samples
  n.samp <- out$mcmc.info$n.samples
  pred <- array(NA, dim=c(500, n.samp,1))
  for(kk in 1:n.samp){
    pred[,kk,1] <- plogis(out$sims.list$mu.alpha[,1][kk] + out$sims.list$mu.alpha[,6][kk] * var.pred.stdz)
  }
  
  #also create vector of occupancy estimate in each station
  site.mean <- mean(psi.covar5)
  site.sd <- sd(psi.covar5)
  site.pred.stdz <- (psi.covar5-site.mean)/site.sd
  
  occu_pred <- plogis(out$mean$mu.alpha[1] + out$mean$mu.alpha[6] * site.pred.stdz)
  
  y.jitter <- jitter(occu_pred, factor=10000) #add noise, so points not overlap
  #the value really depends on the variance of data, so need to calibrate
  
  #then join the biomass prediction with occu_pred
  site_occu <- cbind(psi.covar5, occu_pred)
  names(site_occu)
  colnames(site_occu) <- c("Covariate", "Psi")
  site_occudf <- as.data.frame(site_occu)

  
  #create the plot!
  pmC <- apply(pred, c(1,3), mean)
  criC <- apply(pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
  
  # Plot posterior mean and a random sample of 100 from posterior of regression
  selection <- sample(1:n.samp,500)
  #selection <- sample(1:n.samp,10) #exploratory
  jpeg(file=paste("~", #Directory path
                  "Community Plot Biomass 1000m", ".jpeg", sep=""), width=1000, height=800)
  par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,1.1))
  matplot(var.pred, pred[,selection,], ylab="Mean Occupancy", xlab="Biomass 1000 m", 
          type = "l", lty = 1, lwd=1, col="grey", ylim=c(0,1), xlim=c(covar.min.value, covar.max.value),
          cex.axis=1.5, cex.lab=2)
  lines(var.pred, apply(pred[,,],1,mean), lwd=3, col="blue")
  matlines(var.pred, t(criC[,,]), col="blue", lty=3, lwd=2)
  points(x = site_occudf$Covariate, y= y.jitter, pch = 16,
         col = "black")
  mtext(paste("Community-level Mammalian Occupancy: ", "Biomass", sep=" "), side=3, line=-2.5, cex=1.5, outer=TRUE)
  dev.off()
 }

### Proportion ####
{
  # Create objects to define the parameters of the model prediction
  covar.min.value <- round_any(min(psi.covar6), 1, floor)
  covar.max.value <- round_any(max(psi.covar6), 1, ceiling)
  
  # Calculate covariate values for prediction
  var.pred <- seq(covar.min.value, covar.max.value,,500)
  var.mean <- mean(psi.covar6)
  var.sd <- sd(psi.covar6)
  var.pred.stdz <- (var.pred - var.mean)/var.sd
  
  # Create a temp folder storing species-specific posterior distributions of all monitored parameters
  # Compute predictions based on the number of MCMC samples
  n.samp <- out$mcmc.info$n.samples
  pred <- array(NA, dim=c(500, n.samp,1))
  for(kk in 1:n.samp){
    pred[,kk,1] <- plogis(out$sims.list$mu.alpha[,1][kk] + out$sims.list$mu.alpha[,7][kk] * var.pred.stdz)
  }
  
  #also create vector of occupancy estimate in each station
  site.mean <- mean(psi.covar6)
  site.sd <- sd(psi.covar6)
  site.pred.stdz <- (psi.covar6-site.mean)/site.sd
  
  occu_pred <- plogis(out$mean$mu.alpha[1] + out$mean$mu.alpha[7] * site.pred.stdz)
  
  y.jitter <- jitter(occu_pred, factor=100000) #add noise, so points not overlap
  #the value is so huge to accomodate lack of variance
  
  #then join the biomass prediction with occu_pred
  site_occu <- cbind(psi.covar6, occu_pred)
  names(site_occu)
  colnames(site_occu) <- c("Covariate", "Psi")
  site_occudf <- as.data.frame(site_occu)
  
  
  #create the plot!
  pmC <- apply(pred, c(1,3), mean)
  criC <- apply(pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
  
  # Plot posterior mean and a random sample of 100 from posterior of regression
  selection <- sample(1:n.samp,500)
  #selection <- sample(1:n.samp,10) #exploratory
  jpeg(file=paste("~", #Directory path 
                  "Community Plot Proportion", ".jpeg", sep=""), width=1000, height=800)
  par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,1.1))
  matplot(var.pred, pred[,selection,], ylab="Mean Occupancy", xlab="Proportion", 
          type = "l", lty = 1, lwd=1, col="grey", ylim=c(0,1), xlim=c(covar.min.value, covar.max.value),
          cex.axis=1.5, cex.lab=2)
  lines(var.pred, apply(pred[,,],1,mean), lwd=3, col="blue")
  matlines(var.pred, t(criC[,,]), col="blue", lty=3, lwd=2)
  points(x = site_occudf$Covariate, y=y.jitter,pch = 16,
         col = "black")
  mtext(paste("Community-level Mammalian Occupancy: ", "Proportion", sep=" "), side=3, line=-2.5, cex=1.5, outer=TRUE)
  dev.off()
}

### Distance to water ####
{
  # Create objects to define the parameters of the model prediction
  covar.min.value <- round_any(min(psi.covar7), 1, floor)
  covar.max.value <- round_any(max(psi.covar7), 1, ceiling)
  
  # Calculate covariate values for prediction
  var.pred <- seq(covar.min.value, covar.max.value,,500)
  var.mean <- mean(psi.covar7)
  var.sd <- sd(psi.covar7)
  var.pred.stdz <- (var.pred - var.mean)/var.sd
  
  # Create a temp folder storing species-specific posterior distributions of all monitored parameters
  # Compute predictions based on the number of MCMC samples
  n.samp <- out$mcmc.info$n.samples
  pred <- array(NA, dim=c(500, n.samp,1))
  for(kk in 1:n.samp){
    pred[,kk,1] <- plogis(out$sims.list$mu.alpha[,1][kk] + out$sims.list$mu.alpha[,8][kk] * var.pred.stdz)
  }
  
  #also create vector of occupancy estimate in each station
  site.mean <- mean(psi.covar7)
  site.sd <- sd(psi.covar7)
  site.pred.stdz <- (psi.covar7-site.mean)/site.sd
  
  occu_pred <- plogis(out$mean$mu.alpha[1] + out$mean$mu.alpha[8] * site.pred.stdz)
  
  y.jitter <- jitter(occu_pred, factor=1000) #add noise, so points not overlap
  #the value really depends on the variance of data, so need to calibrate
  
  #then join the biomass prediction with occu_pred
  site_occu <- cbind(psi.covar7, occu_pred)
  names(site_occu)
  colnames(site_occu) <- c("Covariate", "Psi")
  site_occudf <- as.data.frame(site_occu)
  
  
  #create the plot!
  pmC <- apply(pred, c(1,3), mean)
  criC <- apply(pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
  
  # Plot posterior mean and a random sample of 100 from posterior of regression
  selection <- sample(1:n.samp,500)
  #selection <- sample(1:n.samp,10) #exploratory
  jpeg(file=paste("~", #Directory path
                  "Community Plot Water", ".jpeg", sep=""), width=1000, height=800)
  par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,1.1))
  matplot(var.pred, pred[,selection,], ylab="Mean Occupancy", xlab="Distance to water", 
          type = "l", lty = 1, lwd=1, col="grey", ylim=c(0,1), xlim=c(covar.min.value, covar.max.value),
          cex.axis=1.5, cex.lab=2)
  lines(var.pred, apply(pred[,,],1,mean), lwd=3, col="blue")
  matlines(var.pred, t(criC[,,]), col="blue", lty=3, lwd=2)
  points(x = site_occudf$Covariate, y=y.jitter,pch = 16,
         col = "black")
  mtext(paste("Community-level Mammalian Occupancy: ", "Water", sep=" "), side=3, line=-2.5, cex=1.5, outer=TRUE)
  dev.off()
}

### Access ###############
{
  # Create objects to define the parameters of the model prediction
  covar.min.value <- round_any(min(psi.covar8), 1, floor)
  covar.max.value <- round_any(max(psi.covar8), 1, ceiling)
  
  # Calculate covariate values for prediction
  var.pred <- seq(covar.min.value, covar.max.value,,500)
  var.mean <- mean(psi.covar8)
  var.sd <- sd(psi.covar8)
  var.pred.stdz <- (var.pred - var.mean)/var.sd
  
  # Create a temp folder storing species-specific posterior distributions of all monitored parameters
  # Compute predictions based on the number of MCMC samples
  n.samp <- out$mcmc.info$n.samples
  pred <- array(NA, dim=c(500, n.samp,1))
  for(kk in 1:n.samp){
    pred[,kk,1] <- plogis(out$sims.list$mu.alpha[,1][kk] + out$sims.list$mu.alpha[,9][kk] * var.pred.stdz
                          + out$sims.list$mu.alpha[,10][kk] * var.pred.stdz^2)
  }
  
  #also create vector of occupancy estimate in each station
  site.pred.stdz <- (psi.covar8-var.mean)/var.sd
  
  occu_pred <- plogis(out$mean$mu.alpha[1] + out$mean$mu.alpha[9] * site.pred.stdz 
                      + out$mean$mu.alpha[10] * site.pred.stdz^2)
  
  y.jitter <- jitter(occu_pred, factor=1000) #add noise, so points not overlap
  #the value really depends on the variance of data, so need to calibrate
  
  #then join the  prediction with occu_pred
  site_occu <- cbind(psi.covar8, occu_pred)
  names(site_occu)
  colnames(site_occu) <- c("Covariate", "Psi")
  site_occudf <- as.data.frame(site_occu)
  
  
  #create the plot!
  pmC <- apply(pred, c(1,3), mean)
  criC <- apply(pred, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
  
  # Plot posterior mean and a random sample of 100 from posterior of regression
  selection <- sample(1:n.samp,500)
  #selection <- sample(1:n.samp,10) #exploratory
  jpeg(file=paste("~", #Directory path
                  "Community Plot Access", ".jpeg", sep=""), width=1000, height=800)
  par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,1.1))
  matplot(var.pred, pred[,selection,], ylab="Mean Occupancy", xlab="Accessibility", 
          type = "l", lty = 1, lwd=1, col="grey", ylim=c(0,1), xlim=c(covar.min.value, covar.max.value),
          cex.axis=1.5, cex.lab=2)
  lines(var.pred, apply(pred[,,],1,mean), lwd=3, col="blue")
  matlines(var.pred, t(criC[,,]), col="blue", lty=3, lwd=2)
  points(x = site_occudf$Covariate, y=y.jitter,pch = 16,
         col = "black")
  mtext(paste("Community-level Mammalian Occupancy: ", "Accessibility", sep=" "), side=3, line=-2.5, cex=1.5, outer=TRUE)
  dev.off()
}
