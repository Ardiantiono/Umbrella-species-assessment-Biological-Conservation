### Functions to generate occupancy graphs
### Code developed by Nicolas J. Deere adapted by Ardiantiono

rm(list=ls())

# Set working directory (where the all data required to run the models are stored)
setwd("~") #fill with directory path
getwd()

jags.out <- load("~//HMSOM_AllSpecies_DA.RData") #Call the model output 
names(out) #it reads as out not jags.out..

# Create unique identifiers for species, covariate names (shorthand and full), scale and starting covariates
s.list<- c("SP1", "SP2", "SP3", "SP4", "SP5")  

# Create a character string for response variable names
var.names <- c("Terrain Ruggedness Index", "Terrain Ruggedness Index Squared", "Elevation", "Elevation Squared",
               "Aboveground Biomass", "Proportion of Forest Cover",
               "Distance to Water", "Accessibility", "Accessibility Squared")
i.var <- c("TRI", "TRI2", "Elev", "Elev2", "Biomass", "Prop", "d_water", "ttcsm", "ttcsm2") 
cov.i <- 1

var.names.1 <- c("Intercept", "Terrain Ruggedness Index", "Terrain Ruggedness Index Squared", "Elevation", "Elevation Squared",
                "Aboveground Biomass", "Proportion of Forest Cover",
                "Distance to Water", "Accessibility", "Accessibility Squared")

x <- c(2:10) #covariate coefficient from output to be used

##########################################################
#Testing 50% BCI (use this for global model, directory different)

for(i in 1:length(x)) {
  # Get posterior means and 95% BCI for the community and all species
  commD<-out$sims.list$mu.alpha[,x[i]]
  psiD<-out$sims.list$alpha[,,x[i]] #keep only covariate of interest (discard 2nd order polynomial)
  # Calculate mean and 95% credible intervals of posterior distribution
  comm.summ <- quantile(commD, probs=c(0.025,0.25,0.5,0.75,0.975))
  psi.summ <- apply(psiD,2,quantile,prob=c(0.025,0.25,0.5,0.75,0.975)) #2=column, 1=row
  
  # Create graphs looping across all covariates within scale
  # Code colours significant points automatically
  mypath=file.path("~", #Directory path
                   paste("Final catterpilar plot ", i.var[i], ".jpeg", sep=""))
  jpeg(file=mypath, width=600, height=1000)
  par(mfrow=c(1,1),mar=c(2,12,2,2))
  plot(1,1,type="n",xlim=c(-2,2),ylim=c(1,23),xlab="",ylab="",yaxt="n",xaxt="n",cex.main=1.5,cex.lab=1.25)
  axis(1, at = seq(-2, 2, by = 0.5))
  abline(v=comm.summ[1],col="red",lty=2)
  abline(v=comm.summ[5],col="red",lty=2) #can add abline for 25% and 75% if needed
  abline(v=comm.summ[3],lw=3,col="red")
  axis(2,at=1:23,labels=rev(s.list),las=2)
  
  segments(rev(psi.summ[1,]),1:23,rev(psi.summ[5,]),1:23, col="grey")
  arrows(rev(psi.summ[1,]),1:23, rev(psi.summ[5,]),1:23, length=0.05, angle=90, col="grey")
  arrows(rev(psi.summ[5,]),1:23, rev(psi.summ[1,]),1:23, length=0.05, angle=90, col="grey")
  segments(rev(psi.summ[2,]),1:23,rev(psi.summ[4,]),1:23,lw=2,col="grey")
  
  sig1 <- (rev(psi.summ[1,])*rev(psi.summ[5,])) > 0
  segments(rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1],rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], col="blue")
  arrows(rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1], rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], length=0.05, angle=90, col="blue")
  arrows(rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1], length=0.05, angle=90, col="blue")
  segments(rev(psi.summ[2,])[sig1==1],(1:23)[sig1==1],rev(psi.summ[4,])[sig1==1],(1:23)[sig1==1],lw=2,col="blue")
  
  sig2 <- (rev(psi.summ[2,])*rev(psi.summ[4,])) > 0
  segments(rev(psi.summ[2,])[sig2==1],(1:23)[sig2==1],rev(psi.summ[4,])[sig2==1],(1:23)[sig2==1], lwd=2, col="green")
  arrows(rev(psi.summ[2,])[sig2==1],(1:23)[sig2==1], rev(psi.summ[4,])[sig2==1],(1:23)[sig2==1], length=0.05, angle=90, col="green")
  arrows(rev(psi.summ[4,])[sig2==1],(1:23)[sig2==1], rev(psi.summ[2,])[sig2==1],(1:23)[sig2==1], length=0.05, angle=90, col="green")
  
  points(rev(psi.summ[3,]),1:23,pch=19, col="grey")
  points(rev(psi.summ[3,])[sig1==1],(1:23)[sig1==1],pch=19,col="blue")
  
  abline(v=0,lw=1,lty=1,col="black")
  mtext(side=3, line=-1.5, text=paste(var.names[i], sep=" "), cex=1.5, outer=TRUE)
  dev.off()
}


#Caterpillar plot with 95% BCI only
for(i in 1:length(x)) {
    # Get posterior means and 95% BCI for the community and all species
    commD<-out$sims.list$mu.alpha[,x[i]]
    psiD<-out$sims.list$alpha[,,x[i]]  #keep only covariate of interest (discard 2nd order polynomial)
    # Calculate mean and 95% credible intervals of posterior distribution
    comm.summ <- quantile(commD, probs=c(0.025,0.25,0.5,0.75,0.975))
    psi.summ <- apply(psiD,2,quantile,prob=c(0.025,0.25,0.5,0.75,0.975)) #2=column, 1=row
    
    # Create graphs looping across all covariates within scale
    # Code colours significant points automatically
    mypath=file.path("~", #Directory path
                     paste("Catplot test",i.var[i], ".jpeg", sep=""))
    jpeg(file=mypath, width=600, height=1000)
    par(mfrow=c(1,1),mar=c(2,12,2,2))
    plot(1,1,type="n",xlim=c(-2,2),ylim=c(1,23),xlab="",ylab="",yaxt="n",xaxt="n",cex.main=1.5,cex.lab=1.25)
    axis(1, at = seq(-2, 2, by = 0.5))
    abline(v=comm.summ[1],col="red",lty=2)
    abline(v=comm.summ[5],col="red",lty=2)
    abline(v=comm.summ[3],lw=3,col="red")
    axis(2,at=1:23,labels=rev(s.list),las=2)
    segments(rev(psi.summ[1,]),1:23,rev(psi.summ[5,]),1:23, col="grey")
    arrows(rev(psi.summ[1,]),1:23, rev(psi.summ[5,]),1:23, length=0.05, angle=90, col="grey")
    arrows(rev(psi.summ[5,]),1:23, rev(psi.summ[1,]),1:23, length=0.05, angle=90, col="grey")
    segments(rev(psi.summ[2,]),1:23,rev(psi.summ[4,]),1:23,lw=2,col="grey")
    points(rev(psi.summ[3,]),1:23,pch=19, col="grey")
    sig1 <- (rev(psi.summ[1,])*rev(psi.summ[5,])) > 0
    segments(rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1],rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], col="blue")
    arrows(rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1], rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], length=0.05, angle=90, col="blue")
    arrows(rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1], length=0.05, angle=90, col="blue")
    segments(rev(psi.summ[2,])[sig1==1],(1:23)[sig1==1],rev(psi.summ[4,])[sig1==1],(1:23)[sig1==1],lw=2,col="blue")
    points(rev(psi.summ[3,])[sig1==1],(1:23)[sig1==1],pch=19,col="blue")
    abline(v=0,lw=1,lty=1,col="black")
    mtext(side=3, line=-1.5, text=paste(var.names[i], sep=" "), cex=1.5, outer=TRUE)
    dev.off()
  } 

##########################################################
#Plot covariate influences
#########################################################

{
  # Get posterior means and 95% BCI for the community and all species
  commD<-out$sims.list$mu.alpha[,2:10]
  # Calculate mean and 95% credible intervals of posterior distribution
  comm.summ <- apply(commD,2,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
  
  # Create graphs looping across all covariates within scale
  # Code colours significant points automatically
  mypath=file.path("~", #Directory path
                   paste("Final catterpilar plot ", "Covariates", ".jpeg", sep=""))
  jpeg(file=mypath, width=600, height=400)
  par(mfrow=c(1,1),mar=c(2,12,2,2))
  plot(1,1,type="n",xlim=c(-1,1),ylim=c(1,9),xlab="",ylab="",yaxt="n",xaxt="n",cex.main=1.5,cex.lab=1.25)
  axis(1, at = seq(-1, 1, by = 0.5))
  
  axis(2,at=1:9,labels=rev(var.names),las=2)
  
  segments(rev(comm.summ[1,]),1:9,rev(comm.summ[5,]),1:9, col="grey")
  arrows(rev(comm.summ[1,]),1:9, rev(comm.summ[5,]),1:9, length=0.05, angle=90, col="grey")
  arrows(rev(comm.summ[5,]),1:9, rev(comm.summ[1,]),1:9, length=0.05, angle=90, col="grey")
  segments(rev(comm.summ[2,]),1:9,rev(comm.summ[4,]),1:9,lw=2,col="grey")
  
  sig1 <- (rev(comm.summ[1,])*rev(comm.summ[5,])) > 0
  segments(rev(comm.summ[1,])[sig1==1],(1:9)[sig1==1],rev(comm.summ[5,])[sig1==1],(1:9)[sig1==1], col="blue")
  arrows(rev(comm.summ[1,])[sig1==1],(1:9)[sig1==1], rev(comm.summ[5,])[sig1==1],(1:9)[sig1==1], length=0.05, angle=90, col="blue")
  arrows(rev(comm.summ[5,])[sig1==1],(1:9)[sig1==1], rev(comm.summ[1,])[sig1==1],(1:9)[sig1==1], length=0.05, angle=90, col="blue")
  segments(rev(comm.summ[2,])[sig1==1],(1:9)[sig1==1],rev(comm.summ[4,])[sig1==1],(1:9)[sig1==1],lw=2,col="blue")
  
  sig2 <- (rev(comm.summ[2,])*rev(comm.summ[4,])) > 0
  segments(rev(comm.summ[2,])[sig2==1],(1:9)[sig2==1],rev(comm.summ[4,])[sig2==1],(1:9)[sig2==1], lwd=2, col="green")
  arrows(rev(comm.summ[2,])[sig2==1],(1:9)[sig2==1], rev(comm.summ[4,])[sig2==1],(1:9)[sig2==1], length=0.05, angle=90, col="green")
  arrows(rev(comm.summ[4,])[sig2==1],(1:9)[sig2==1], rev(comm.summ[2,])[sig2==1],(1:9)[sig2==1], length=0.05, angle=90, col="green")
  
  points(rev(comm.summ[3,]),1:9,pch=19, col="grey")
  points(rev(comm.summ[3,])[sig1==1],(1:9)[sig1==1],pch=19,col="blue")
  
  abline(v=0,lw=1,lty=1,col="black")
  mtext(side=3, line=-1.5, text=paste("Covariates Influence", sep=" "), cex=1.5, outer=TRUE)
  dev.off()
}

###For elevation only, if needed to increase axis scale 

{
  # Get posterior means and 95% BCI for the community and all species
  commD<-out$sims.list$mu.alpha[,x[3]] #no 3, not 4 I wonder why...
  psiD<-out$sims.list$alpha[,,x[3]] #keep only covariate of interest (discard 2nd order polynomial)
  # Calculate mean and 95% credible intervals of posterior distribution
  comm.summ <- quantile(commD, probs=c(0.025,0.25,0.5,0.75,0.975))
  psi.summ <- apply(psiD,2,quantile,prob=c(0.025,0.25,0.5,0.75,0.975)) #2=column, 1=row
  
  # Create graphs looping across all covariates within scale
  # Code colours significant points automatically
  mypath=file.path("~", #Directory path
                   paste("Final catterpilar plot Elevation", ".jpeg", sep=""))
  jpeg(file=mypath, width=800, height=1000)
  par(mfrow=c(1,1),mar=c(2,12,2,2))
  plot(1,1,type="n",xlim=c(-6,4),ylim=c(1,23),xlab="",ylab="",yaxt="n",xaxt="n",cex.main=1.5,cex.lab=1.25)
  axis(1, at = seq(-6, 4, by = 1))
  abline(v=comm.summ[1],col="red",lty=2)
  abline(v=comm.summ[5],col="red",lty=2) #can add abline for 25% and 75% if needed
  abline(v=comm.summ[3],lw=3,col="red")
  axis(2,at=1:23,labels=rev(s.list),las=2)
  
  segments(rev(psi.summ[1,]),1:23,rev(psi.summ[5,]),1:23, col="grey")
  arrows(rev(psi.summ[1,]),1:23, rev(psi.summ[5,]),1:23, length=0.05, angle=90, col="grey")
  arrows(rev(psi.summ[5,]),1:23, rev(psi.summ[1,]),1:23, length=0.05, angle=90, col="grey")
  segments(rev(psi.summ[2,]),1:23,rev(psi.summ[4,]),1:23,lw=2,col="grey")
  
  sig1 <- (rev(psi.summ[1,])*rev(psi.summ[5,])) > 0
  segments(rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1],rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], col="blue")
  arrows(rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1], rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], length=0.05, angle=90, col="blue")
  arrows(rev(psi.summ[5,])[sig1==1],(1:23)[sig1==1], rev(psi.summ[1,])[sig1==1],(1:23)[sig1==1], length=0.05, angle=90, col="blue")
  segments(rev(psi.summ[2,])[sig1==1],(1:23)[sig1==1],rev(psi.summ[4,])[sig1==1],(1:23)[sig1==1],lw=2,col="blue")
  
  sig2 <- (rev(psi.summ[2,])*rev(psi.summ[4,])) > 0
  segments(rev(psi.summ[2,])[sig2==1],(1:23)[sig2==1],rev(psi.summ[4,])[sig2==1],(1:23)[sig2==1], lwd=2, col="green")
  arrows(rev(psi.summ[2,])[sig2==1],(1:23)[sig2==1], rev(psi.summ[4,])[sig2==1],(1:23)[sig2==1], length=0.05, angle=90, col="green")
  arrows(rev(psi.summ[4,])[sig2==1],(1:23)[sig2==1], rev(psi.summ[2,])[sig2==1],(1:23)[sig2==1], length=0.05, angle=90, col="green")
  
  points(rev(psi.summ[3,]),1:23,pch=19, col="grey")
  points(rev(psi.summ[3,])[sig1==1],(1:23)[sig1==1],pch=19,col="blue")
  
  abline(v=0,lw=1,lty=1,col="black")
  mtext(side=3, line=-1.5, text=paste("Elevation", sep=" "), cex=1.5, outer=TRUE)
  dev.off()
}
