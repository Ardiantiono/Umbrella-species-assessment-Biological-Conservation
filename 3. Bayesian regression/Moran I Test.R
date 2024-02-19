# Test Moran I Test for biodiversity parameters
# Source https://stats.oarc.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/

# Clear workspace
rm(list=ls())

# Load required library
library(ape)

# Set working directory (where the all data required to run the models are stored)
setwd("D:\\Data analysis\\Chapter 1_Surrogate analysis\\6. Bayesian regression")
getwd()

# Load data
data <- read.csv("Species_Biodiversity_Leuser - SES Z.csv", header=TRUE)
head(data)


# Generate inverse distance matrix using lat long coordinates
data.dists <- as.matrix(dist(cbind(data$Longitude, data$Latitude)))
data.dists <- data.dists/1000 #scale to km

data.dists.inv <- 1/data.dists
diag(data.dists.inv) <- 0
summary(data.dists.inv)



# Calculate Moran's I test for biodiversity parameter

Moran.I(data$Occu_Mean, data.dists.inv)
Moran.I(data$SPR_Mean, data.dists.inv)
Moran.I(data$SESFDis, data.dists.inv)
Moran.I(data$mpd.obs.z, data.dists.inv)
Moran.I(data$mntd.obs.z, data.dists.inv)

library(ggplot2) #clear clustering, SPR less obvious
p <- ggplot(data, aes(Longitude, Latitude))
p + geom_point(aes(colour = Occu_Mean))

p <- ggplot(data, aes(Longitude, Latitude))
p + geom_point(aes(colour = SPR_Mean))

p <- ggplot(data, aes(Longitude, Latitude))
p + geom_point(aes(colour = SESFDis))

p <- ggplot(data, aes(Longitude, Latitude))
p + geom_point(aes(colour = mpd.obs.z))

p <- ggplot(data, aes(Longitude, Latitude))
p + geom_point(aes(colour = mntd.obs.z))


### Plot spline colleogram
library(ncf)

occuI <- spline.correlog(x= data$Longitude,
                         y= data$Latitude,
                         z= data$Occu_Mean,
                         resamp= 100,
                         xmax= 40000,
                         quiet= TRUE)

plot(occuI)

SRI <- spline.correlog(x= data$Longitude,
                         y= data$Latitude,
                         z= data$SPR_Mean,
                         resamp= 100,
                         xmax= 40000,
                         quiet= TRUE)

plot(SRI)

FDI <- spline.correlog(x= data$Longitude,
                       y= data$Latitude,
                       z= data$SESFDis,
                       resamp= 100,
                       xmax= 40000,
                       quiet= TRUE)

plot(FDI)

MPDI <- spline.correlog(x= data$Longitude,
                       y= data$Latitude,
                       z= data$mpd.obs.z,
                       resamp= 100,
                       xmax= 40000,
                       quiet= TRUE)

plot(MPDI)

MNTDI <- spline.correlog(x= data$Longitude,
                       y= data$Latitude,
                       z= data$mntd.obs.z,
                       resamp= 100,
                       xmax= 40000,
                       quiet= TRUE)

plot(MNTDI)

#Arrange multiple plot
library(gridExtra)
library(grid)


#Save as picture
png(file="Spline colleogram Moran's I Test.png",
    width=600, height=500)
par(mfrow = c(2, 2), mai = c(0.4, 0.4, 0.3, 0.1)) #bottom, left, top, right
plot(occuI)
plot(SRI)
plot(FDI)
plot(MPDI)
dev.off()

par(mfrow = c(1, 1))
