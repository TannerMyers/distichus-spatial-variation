# Custom function from https://github.com/LauraANunes/NunesPearson2017
# From: Nunes, L. A. and Pearson, R. G. (2017), 
# A null biogeographical test for assessing ecological niche evolution.
# Journal of Biogeography. doi: 10.1111/jbi.12910

### Trimming occurences, remove occurences which represent outliers in environmental hyperspace (Farber & Kadmon, 2003)
trim.occ<-function(bios,occ, plot=TRUE){
  S=matrix(c(extract(bios,occ)),ncol=nlayers(bios),nrow=nrow(occ))
  C=cov(S) #covariance matrix
  M=apply(S,2,mean) #mean conditions of each climate indice
  
  D<-rep(NA,nrow(S))
  for(i in 1:nrow(S)){
    
    Tr=t(S[i,]-M) #transpose opeerator
    D[i]=Tr%*%((C)^-1)%*%(S[i,]-M) #malahanobis distance
  }
  
  D<-cbind(c(1:nrow(S)),D) #add locality index to D distances
  
  qD<-D[c(which(D[,2]>quantile(D,0.05))),] #trim 5%
  qD<-qD[c(which(qD[,2]<quantile(D,0.95))),] #trim 95%
  
  occ_trim<-occ[qD[,1],] #keep locality points within 5-95th percentile
  if(plot==TRUE){
    plot(D[,2], xlab='', ylab='Mahalanobis Distance')
    abline(h=quantile(D,0.95), col='red')
    abline(h=quantile(D,0.05), col='red')
    plot(bios[[1]], col='grey')
    points(occ)
    points(occ_trim,col='red')
  }
  return(occ_trim)
}