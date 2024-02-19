## Calculate MPD and MNTD and their Standardized Effect Size ##

#Getting Phylogenetic Tree for Mammals in Leuser ##
# Set working directory (where the all data required to run the models are stored)
rm(list=ls())
setwd("~") #Set up file path
getwd()


#########################################
library(ape)

# Loading Mammalian Supertree
mammal.tree <- read.nexus(file="ELE_1307_sm_SA1.txt")
mammal.tree.1 <- mammal.tree[[1]] 

class(mammal.tree)

# Pruning mammalian supertree
sp <- read.csv("Species_list.csv", header = FALSE)
sp_list<-sp[,1]

library(geiger)
dat <- matrix(data=NA, nrow=length(sp_list), ncol=1,
              dimnames=list(sp_list, NA))
head(dat)
View(dat)
tree.leuser <- treedata(mammal.tree.1, dat)$phy
tree.leuser 
class(tree.leuser)


# Plotting phylogenetic tree
plot.phylo(tree.leuser, cex=1, no.margin=TRUE) #If worked should show the phylo tree
help(plot.phylo)


#####################################
## Calculate MPD and MNTD + Standardized Effect Size
## See https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf as reference

library(picante)

phydist <- cophenetic(tree.leuser) #convert to distance matrix

## Calculate MPD using "richness" algorithm
ses.mpd.result <- ses.mpd(community, phydist, null.model="richness",
                          abundance.weighted=TRUE, runs=1000)
ses.mpd.result

# Result: ntaxa (richness), mpd.obs (mpd value), mpd.obs.z (SES), mpd.obs.p (p-value)


## Calculate MNTD
ses.mntd.result <- ses.mntd(community, phydist, null.model="richness", 
                              abundance.weighted=TRUE, runs=1000)
ses.mntd.result


#Save as file
#Combine datasets
SES.PD <- cbind.data.frame(ses.mpd.result, ses.mntd.result)
write.csv(SES.PD, file=paste("XXXX", ".csv", sep="")) #assign the file name



