# Algorithm for determining how to allocate additional samples within a study area given some existing samples
# The method is based mainly on the Carre et al (2007) HELS and  algorithm
# Hypercube Evaluation of a Legacy Sample (HELS)
#
# basically it entails:
# 1. From the covariate data generate a matrix of quantiles (n x p) n = number of samples that are needed. p = number of available covariates
# 2. Using the quantile matrix create Hypercube matrices for both the existing point data and covariates.
# 3. Work out the densities of the elements in the hypercubes (count of obs/pixels within each quantile / total number of data or pixels)
# 4. Work out the ratio of sampling and grid densities... This will identify under and over sampling in the hypercube
# 5. To add additional samples:
#         1. rank the ratios from smallest to largest
#         2. workout the number of samples required to equalise quantile density of grids and sample data
#         3. Repeat step 5.2 until total number of additonal samples have been allocated.
## Created: 18/05/17

setwd("Z:/Dropbox/2019/rmuddles/clhs_addtion")

# Libraries
library(raster);library(sp); library(rgdal)


# Number of addtional samples to take
nosP<- 100



##INPUT DATA
#########################################################################################################################
#covariates
tempD<- readRDS("tempD.rds")
str(tempD)
names(tempD)

# rasterise data
s1<- stack()
for (i in 4:ncol(tempD)){
  r1<- rasterFromXYZ(tempD[,c(1,2,i)])
  names(r1)<- names(tempD[i])
  s1<- stack(s1,r1)
}
s1



#quantile matrix (of the covariate data)
# number of quantiles = nosP
# note there are 7 covariate layers. First covariate layer in column 4
q.mat<- matrix(NA, nrow=(nosP+1), ncol= 7)
j=1
for (i in 4:ncol(tempD)){ #not the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(tempD[,i]) - min(tempD[,i])
  step1<- ran1/nosP 
  q.mat[,j]<- seq(min(tempD[,i]), to = max(tempD[,i]), by =step1)
  j<- j+1}
#q.mat




#count of pixels within each quantile for each covariate
#############################################
## This takes a while to do so only do it once if you can 

cov.mat<- matrix(0, nrow=nosP, ncol=7)

for (i in 1:nrow(tempD)){ # the number of pixels
  cntj<- 1 
  for (j in 4:ncol(tempD)){ #for each column
    dd<- tempD[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}

#saveRDS(cov.mat, file = "ancil_cov.rds")
#readRDS("ancil_cov.rds")

cov.mat[which(cov.mat==0)]<- 0.000001 # very small number so we dont have to deal with zeros
# covariate data density
covDens<- cov.mat/nrow(tempD) 
#covDens # covariate data density matrix
################################################################################



################################################################################
#Point data (the exisitng points)
dat<- read.table("HunterValley_SiteObsAll.txt", header = T,sep = ",")  # existing soil point data


#extract covariate data at points
coordinates(dat)<- ~ X +Y 
DSM_data<- extract(s1,dat, sp= 1, method = "simple") #extract
dat<- as.data.frame(DSM_data)
dat<- dat[complete.cases(dat),]

#Point dat (little data manipulations)
dat<- dat[,c(2,3,4,6:ncol(dat))]
names(dat)[1:3]<- c("x", "y", "cellNos")
dat$cellNos<- 0
str(dat)


# Count of sample sample data for each covariate quantile
# essentially the same as what was done for the gridded covariate data
h.mat<- matrix(0, nrow=nosP, ncol=7)

for (i in 1:nrow(dat)){ # the number of observations
  cntj<- 1 
  for (j in 4:ncol(dat)){ #for each column
    dd<- dat[i,j]  
    for (k in 1:nosP){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[k, cntj]<- h.mat[k, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

#density
h.mat[which(h.mat==0)]<- 0.000001  # small number so we dont have to deal with zeros
datDens<-h.mat/nrow(dat) # data density


###########################################################################################3



############################################################################################
#### selecting new samples

rat<- datDens/covDens # ratio of data density and covariate density
or<- order(rat) # rank the index where the biggest discrepancy is at the top


## indexes of quantiles that are not adequately sampled ie where rat is less than 1
l1<- which(rat < 1, arr.ind = T)
l1<- cbind(l1,which(rat < 1) )
length(rat) # total number of quantiles
nrow(l1) # number of quanitles where there is a discrepancy



#start from the highest discrepancy then work our way down
upSamp<- nosP
rpos<- 1
base<- 3 # constant so that covariate columns are easily to select (related to the column positions of the covariates in the tempD data frame)


while (upSamp != 0){  # while the number of samples to allocate is greater than 0
  
  # quantile position where the discrepancy is
  indy<- which(l1[,3]==or[rpos]) 
  rp<- l1[indy, 1]
  rc<- l1[indy, 2]
  
  ex<- floor(nrow(dat) * (datDens[rp,rc])) #existing count of samples within the selected quantile
  eq<- ceiling(nrow(dat) * (covDens[rp,rc])) # number of samples needed to get to equal density between data and covariates
  sn<- eq-ex #number of samples needed
  if (upSamp < sn) {sn <- upSamp} # just so we dont over allocate
  
  
  #covariate selection
  covL<- q.mat[rp, rc]
  covU<- q.mat[rp+1, rc]
  subDat<- tempD[tempD[,(base+rc)] >= covL & tempD[,(base+rc)] <= covU,] # subset the covariates that meet the standard
  
  training <- sample( nrow(subDat), sn) #random number
  subDat2<- subDat[training,]
  
  #remove selected samples from tempD so that repeated sampling does not occur (Is this necessary??)
  tempD<- tempD[!(tempD$cellNos %in% subDat2$cellNos), ]
  
  # Append new data to sampling dataframe
  dat<- rbind(dat,subDat2)
  
  #adjust the while params
  rpos<- rpos + 1 # move to next largest discrepancy
  upSamp<- upSamp - sn # update the required sample number
  print(sn)}




### Checking sample data against coobs

# Specify the different surveys (original and addtional)
dat$survey<- NA
str(dat)
dat[dat$cellNos==0,"survey"]<- 1
dat[dat$cellNos!=0,"survey"]<- 2

saveRDS(dat, file = "hels_datout.rds")

## Spatial points
coordinates(dat) <- ~x + y
str(dat)

## Coordinate reference systems
proj4string(dat) <- CRS("+init=epsg:32756")
dat@proj4string

## Write point data to shapefile
writeOGR(dat, ".", "HV_dat_shape", "ESRI Shapefile")





# The following raster has been derived by comparing the similarity between the mulivariate values of the grids and observed data
# Low numbers mean that there are not many data points showing similarity to the grid cell location. High mumber means quite a few
# obervations are similar

#raster extraction
#covariates
r1<- raster("sampleNos.tif")
DSM_data2<- extract(r1,dat, sp= 1, method = "simple")



# Doing some summary statistics between the raster grid values and the sample sites for both original and addtional data
#Exisitng sample data
dat1<- DSM_data2[DSM_data2$survey ==1, ]
sum(dat1$sampleNos >= 0 & dat1$sampleNos <= 5) / nrow(dat1)
sum(dat1$sampleNos > 5 & dat1$sampleNos <= 10) / nrow(dat1)
sum(dat1$sampleNos > 10 & dat1$sampleNos <= 20) / nrow(dat1)
sum(dat1$sampleNos > 20 & dat1$sampleNos <= 40) / nrow(dat1)
sum(dat1$sampleNos > 40) / nrow(dat1)

# additional sample
dat2<- DSM_data2[DSM_data2$survey !=1, ]
sum(dat2$sampleNos >= 0 & dat2$sampleNos <= 5) / nrow(dat2)
sum(dat2$sampleNos > 5 & dat2$sampleNos <= 10) / nrow(dat2)
sum(dat2$sampleNos > 10 & dat2$sampleNos <= 20) / nrow(dat2)
sum(dat2$sampleNos > 20 & dat2$sampleNos <= 40) / nrow(dat2)
sum(dat2$sampleNos > 40) / nrow(dat2)

save.image("clhs_samp_hels.RData") #save R session



