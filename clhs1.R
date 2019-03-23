# Conditioned Latin Hypercube sampling
# Checking how the clhs R function handles existing sample data 
# What does 'include' do
# Checking outputs with coobs map to detemine whether new samples go to areas of low environmental coverage
# created: 15.3.2019



setwd("Z:/Dropbox/2019/rmuddles/clhs_addtion")

# data frame of covariate data
tempD<- readRDS("tempD.rds")



# Libraries
library(raster);library(sp); library(rgdal); library(clhs)



#rasterise covariate data
s1<- stack()
for (i in 4:ncol(tempD)){
  r1<- rasterFromXYZ(tempD[,c(1,2,i)])
  names(r1)<- names(tempD[i])
  s1<- stack(s1,r1)
}
s1


#tabulate (as i want the cell number)
tempD <- data.frame(cellNos = seq(1:ncell(s1)))
vals <- as.data.frame(getValues(s1))
tempD<- cbind(tempD, vals)
tempD <- tempD[complete.cases(tempD), ]
cellNos <- c(tempD$cellNos)
gXY <- data.frame(xyFromCell(s1, cellNos, spatial = FALSE))
tempD<- cbind(gXY, tempD)
str(tempD)

#rasterise again
s1<- stack()
for (i in 3:ncol(tempD)){
  r1<- rasterFromXYZ(tempD[,c(1,2,i)])
  names(r1)<- names(tempD[i])
  s1<- stack(s1,r1)
}
s1



#Point data (the exisitng sample point data)
dat<- read.table("HunterValley_SiteObsAll.txt", header = T,sep = ",")  # existing soil point data


#extract covariate data at points
coordinates(dat)<- ~ X +Y 
DSM_data<- raster::extract(s1,dat, sp= 1, method = "simple") #extract
dat<- as.data.frame(DSM_data)
dat<- dat[complete.cases(dat),]
str(dat)


#remove the grid points where there is point data
tempD<- tempD[-which(tempD$cellNos %in% dat$cellNos),]
str(tempD)

## combine grid data with the observed data
str(dat)
dat.sub<- dat[,c(2,3,6:13)]
names(dat.sub)[1:2]<- c("x", "y")
tempD.new<- rbind(dat.sub, tempD)
tempD.new$type<- NA
tempD.new$type[1:nrow(dat.sub)]<- "orig"
tempD.new$type[(nrow(dat.sub)+1):nrow(tempD.new)]<- "possibles"

## clhs sampling with fixed obs and add an extra 100 sites
# note usage of the include parameter
names(tempD.new)
res <- clhs(tempD.new[,c(4:10)], size = nrow(dat.sub) + 100, 
            iter = 10000, progress = TRUE, simple = TRUE, 
            include = c(1:nrow(dat.sub)))
res
saveRDS(res, file = "clhs_res.rds")


# get the selected data
dat.sel<- tempD.new[res,]


## coobs raster surface
r1<- raster("sampleNos.tif")
plot(r1)


# extract coobs dat
coordinates(dat.sel)<- ~ x + y
DSM_data2<- extract(r1,dat.sel, sp= 1, method = "simple")



# Doing some summary statistics between the raster grid values and the sample sites for both original and addtional data
# original sample sites
dat1<- DSM_data2[101:432, ]
sum(dat1$sampleNos >= 0 & dat1$sampleNos <= 5) / nrow(dat1) # very low coobs
sum(dat1$sampleNos > 5 & dat1$sampleNos <= 10) / nrow(dat1) # low coobs
sum(dat1$sampleNos > 10 & dat1$sampleNos <= 20) / nrow(dat1) # moderate coobs
sum(dat1$sampleNos > 20 & dat1$sampleNos <= 40) / nrow(dat1) # high coobs
sum(dat1$sampleNos > 40) / nrow(dat1) # quite high coobs

# additional data as selected by clhs
dat2<- DSM_data2[1:100, ]
sum(dat2$sampleNos >= 0 & dat2$sampleNos <= 5) / nrow(dat2)
sum(dat2$sampleNos > 5 & dat2$sampleNos <= 10) / nrow(dat2)
sum(dat2$sampleNos > 10 & dat2$sampleNos <= 20) / nrow(dat2)
sum(dat2$sampleNos > 20 & dat2$sampleNos <= 40) / nrow(dat2)
sum(dat2$sampleNos > 40) / nrow(dat2)

save.image("clhs_samp.RData") #save R session



