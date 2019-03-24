## clhs work
## KL-divergence calculation for clhs sample with addtional samples and with adapted hels

library(entropy)

setwd("Z:/Dropbox/2019/rmuddles/clhs_addtion")


# data sets
tempD<- readRDS("tempD.rds")

## clhs data
clhs.dat<- readRDS("clhs_datout.rds")


## adapted hels
hels.dat<- readRDS("hels_datout.rds")




# Number of bins
nb<- 50

#quantile matrix (of the covariate data)
str(tempD)
q.mat<- matrix(NA, nrow=(nb+1), ncol= 7)
j=1
for (i in 4:ncol(tempD)){ #note the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(tempD[,i]) - min(tempD[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(tempD[,i]), to = max(tempD[,i]), by =step1)
  j<- j+1}
q.mat



#covariate data hypercube 
## This takes a while to do so only do it once if you can 

cov.mat<- matrix(1, nrow=nb, ncol=7)
for (i in 1:nrow(tempD)){ # the number of pixels
  cntj<- 1 
  for (j in 4:ncol(tempD)){ #for each column
    dd<- tempD[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}

cov.mat



## Hypercube for clhs dat
str(clhs.dat)
h.mat.clhs<- matrix(1, nrow=nb, ncol=7)

for (ii in 1:nrow(clhs.dat)){ # the number of observations
  cntj<- 1 
  for (jj in 4:10){ #for each column
    dd<- clhs.dat[ii,jj]  
    for (kk in 1:nb){  #for each quantile
      kl<- q.mat[kk, cntj] 
      ku<- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat.clhs[kk, cntj]<- h.mat.clhs[kk, cntj] + 1}
    }
    cntj<- cntj+1}}
    
h.mat.clhs
      


## Hypercube for adapted HELS data
str(hels.dat)
h.mat.hels<- matrix(1, nrow=nb, ncol=7)
      
for (ii in 1:nrow(hels.dat)){ # the number of observations
  cntj<- 1 
  for (jj in 4:10){ #for each column
    dd<- hels.dat[ii,jj]  
    for (kk in 1:nb){  #for each quantile
      kl<- q.mat[kk, cntj] 
      ku<- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat.hels[kk, cntj]<- h.mat.hels[kk, cntj] + 1}
    }
    cntj<- cntj+1}}


h.mat.hels





#Kullback-Leibler (KL) divergence
# clhs data
klo.1<- KL.empirical(c(cov.mat[,1]), c(h.mat.clhs[,1])) #1
klo.2<- KL.empirical(c(cov.mat[,2]), c(h.mat.clhs[,2])) #2
klo.3<- KL.empirical(c(cov.mat[,3]), c(h.mat.clhs[,3])) #3
klo.4<- KL.empirical(c(cov.mat[,4]), c(h.mat.clhs[,4])) #4
klo.5<- KL.empirical(c(cov.mat[,5]), c(h.mat.clhs[,5])) #5
klo.6<- KL.empirical(c(cov.mat[,6]), c(h.mat.clhs[,6])) #6
klo.7<- KL.empirical(c(cov.mat[,7]), c(h.mat.clhs[,7])) #7
klo.clhs<- mean(c(klo.1, klo.2,klo.3,klo.4,klo.5,klo.6,klo.7))
klo.clhs

# hels data
klo.1<- KL.empirical(c(cov.mat[,1]), c(h.mat.hels[,1])) #1
klo.2<- KL.empirical(c(cov.mat[,2]), c(h.mat.hels[,2])) #2
klo.3<- KL.empirical(c(cov.mat[,3]), c(h.mat.hels[,3])) #3
klo.4<- KL.empirical(c(cov.mat[,4]), c(h.mat.hels[,4])) #4
klo.5<- KL.empirical(c(cov.mat[,5]), c(h.mat.hels[,5])) #5
klo.6<- KL.empirical(c(cov.mat[,6]), c(h.mat.hels[,6])) #6
klo.7<- KL.empirical(c(cov.mat[,7]), c(h.mat.hels[,7])) #7
klo.hels<- mean(c(klo.1, klo.2,klo.3,klo.4,klo.5,klo.6,klo.7))
klo.hels


