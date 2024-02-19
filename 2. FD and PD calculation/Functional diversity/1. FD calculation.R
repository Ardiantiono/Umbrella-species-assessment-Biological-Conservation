### Calculating FD and Standardized Effect Size ###
### Using true occupancy (Z) instead of psi ###
### Code adapted from Gorczynski et al. (https://royalsocietypublishing.org/doi/10.1098/rspb.2020.2098)
### and Penjor et al  (https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13613)###

rm(list=ls())
# Set working directory (where the all data required to run the models are stored)
setwd("~") #Set up file path
getwd()


#############################################
##Calculate FD (e.g. FRich and FDis) for all sites
##both realized community and species pools
library(tidyr)
library(dplyr)

##Import species traits
spp.pool.traits <- read.csv("Trait_species.csv") 
names(spp.pool.traits)

#Rename rows as species name
rownames(spp.pool.traits) <- spp.pool.traits$Species
spp.pool.traits$Species <- NULL

##Import Species pool list
Species1 <- read.csv("Species_list_FD")

#Import Occupancy data (here we use mean of Z or true occurence, instead of psi as suggested by Jarzyna et al)
#We round the decimal to 2, so very low z would be assigned as zero (or absence)
All.TEAM.Occupancy <- read.csv("Species_True_Occu_FD.csv")
colnames(All.TEAM.Occupancy) <- c("Site","Species","2017")
Occ <- gather(All.TEAM.Occupancy, "Year", "Occupancy", 3)


#weighted by occupancy value 
Species2 <- Occ
Species2$Community <- paste(Species2$Site, Species2$Year) #create new ID
Species2 <- Species2[,-c(1,3)] #retain only species, occupancy, and community
Species2 <- Species2 %>% spread(Community,Occupancy)
Species2[is.na(Species2)] <- 0
Species2 <- filter(Species2, Species2$Species %in% rownames(spp.pool.traits))
rownames(Species2) <- Species2$Species
Species2$Species <- NULL
Species2 <- t(Species2)
Species2 <- Species2[which(rowSums(Species2) > 0),] 


#standardize covariate values
spp.pool.traits <- spp.pool.traits[order(rownames(spp.pool.traits)),]
spp.pool.traits$Bodymass_kg <- log(spp.pool.traits$Bodymass_kg)
spp.pool.traits$Max_longevity_yr <- log(spp.pool.traits$Max_longevity_yr)


#Calculate Functional Diversity
library(FD)
x <- dbFD(spp.pool.traits, Species2, w = c(5,5,5,5,5,5,15,15,15,15)) #weight the trait column
#FRic has different values now as we have 0 values for some species in site-level

fdis <- as.data.frame(x["FDis"])
fric <- as.data.frame(x["FRic"])
feve <- as.data.frame(x["FEve"])
fdiv <- as.data.frame(x["FDiv"])

FD.df <- cbind(fdis, fric, feve, fdiv)

write.csv(FD.df, file=paste("XXX", ".csv", sep="")) #assign the file name


#####################################
## Calculate Standardized Effect Size

library(picante)

# Null model permutations
set.seed(1)

null.model.type <- 'richness' #independentswap works, just small change...
nperm <- 1000

FRic <- matrix(NaN, nrow=nperm, ncol=length(x$FRic))
FDiv <- matrix(NaN, nrow=nperm, ncol=length(x$FDiv))
FDis <- matrix(NaN, nrow=nperm, ncol=length(x$FDis))
FEve <- matrix(NaN, nrow=nperm, ncol=length(x$FEve))

abund.null <- vector(mode='list', nperm)


t1 <- Sys.time()
for(i in  seq(nperm)) { # loop for each permutation
  # Create permuted abundance matrix
  abund.i <- randomizeMatrix(Species2, null.model=null.model.type)
  traits.i <- spp.pool.traits
  
  abund.null[[i]] <- abund.i
  
  # Record results
  res.i <- dbFD(traits.i, abund.i, w = c(5,5,5,5,5,5,15,15,15,15))
  
  # Save results
  FRic[i, ] <- res.i$FRic
  FDiv[i, ] <- res.i$FDiv
  FDis[i, ] <- res.i$FDis
  FEve[i, ] <- res.i$FEve
  
  rm(abund.i); rm(res.i)
  print(paste0(round(i/nperm*100), '% of permutations completed'))
}

t2 <- Sys.time()
t2 - t1  # time elapsed - 7.68 mins for 10 perm
(t2 - t1) / nperm * 1000 # time required for 1000 permutations (10 hours)

# Summary stats
# Standardised effect size
# FDis 
SESFDis = (x$FDis - apply(FDis, 2, mean)) / apply(FDis, 2, sd)  # standardised effect size
qFDis <- NaN * x$FDis
for(i in seq(qFDis)) {
  qFDis[i] <- sum(x$FDis[i] > FDis[,i]) / length(FDis[,i])
}
sigFDis <- qFDis < 0.05 | qFDis > 0.95 # test if outside distribution

SES.FDis <- as.data.frame(SESFDis)
q.FDis <- as.data.frame(qFDis)

# Summary stats
# Standardised effect size
# FRic 
SESFRic = (x$FRic - apply(FRic, 2, mean)) / apply(FRic, 2, sd)  # standardised effect size
qFRic <- NaN * x$FRic
for(i in seq(qFRic)) {
  qFRic[i] <- sum(x$FRic[i] > FRic[,i]) / length(FRic[,i])
}
sigFRic <- qFRic < 0.05 | qFRic > 0.95 # test if outside distribution

SES.FRic <- as.data.frame(SESFRic)
q.FRic <- as.data.frame(qFRic)

#Save as file
#Combine datasets
SES.FD <- cbind.data.frame(SES.FDis, q.FDis, SES.FRic, q.FRic)
write.csv(SES.FD, file=paste("XXX", ".csv", sep="")) #assign the file name


########################################
### Plot the trait distribution (following Gorczynski et al.)

##Figure 1###
##Plots of first three dimmensions of FRich for all sites 
library(tidyverse)
library(FactoMineR)
library(geometry)
library(rgl)
library(FD)

spp.pool.traits <- read.csv("Trait_species.csv", row.names = 1)
spp.pool.traits$Bodymass_kg <- log(spp.pool.traits$Bodymass_kg)
spp.pool.traits$Max_longevity_yr <- log(spp.pool.traits$Max_longevity_yr)
x.dist <- gowdis(spp.pool.traits, w = c(5,5,5,5,5,5,15,15,15,15))
x.dist2 <- sqrt(x.dist)
x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
#Spp.PCA <- x.pco$li[,1:2]
Spp.PCA <- x.pco$li[,1:3]
Spp.PCA <- as.data.frame(Spp.PCA)
Spp.PCA$Species <- rownames(Spp.PCA)

##All species plot##
site.conv <- convhulln(Spp.PCA[1:3], output.options = TRUE)
#rgl.open()
open3d()
#rgl.bg(color = "white")
bg3d(color = "white")
plot(site.conv, color = "grey", alpha = 0.2)
#rgl.points(Spp.PCA[,1],Spp.PCA[,2],Spp.PCA[,3], color = "black", size = 2)
points3d(Spp.PCA[,1],Spp.PCA[,2],Spp.PCA[,3], color = "black", size = 4)


#if want to highlight some species
points3d(0.188384258, 0.255915486, 0.07649959, color = "green", size = 7) #Sp1
points3d(-0.228877728, 0.152305219, -0.05944757, color = "yellow", size = 7) #Sp2
points3d(0.170875071, 0.341858045, 0.02529553, color = "red", size = 7) #Sp3
points3d(-0.460955233, 0.070483146, 0.08035392, color = "blue", size = 7) #S4
points3d(-0.080591988, -0.163488802, -0.24553194, color = "pink", size = 7) #Sp5


segments3d(c(-.4, .4), c(0, 0), c(0, 0), color = "black")
segments3d(c(0, 0), c(-.4,.4), c(0, 0), color = "black")
segments3d(c(0, 0), c(0, 0), c(-.4,.4), color = "black")
#rgl.postscript("All_Species.pdf",fmt="pdf")
rgl.snapshot("All_Species.png", fmt ="png")


