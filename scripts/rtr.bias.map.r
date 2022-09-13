# Custom function from https://github.com/LauraANunes/NunesPearson2017
# From: Nunes, L. A. and Pearson, R. G. (2017), 
# A null biogeographical test for assessing ecological niche evolution.
# Journal of Biogeography. doi: 10.1111/jbi.12910

#### function to measure and visualise spatial bias ###
rtr.bias.map<-function(sp1,sp2,bios,iter,plot=TRUE){
  
  r<-raster(ncol=ncol(bios[[1]]),nrow=nrow(bios[[1]]))
  extent(r)<-extent(bios[[1]])
  res(r)<-res(bios[[1]])
  projection(r)<-projection(bios[[1]])
  values(r)<-0
  r[which(is.na(values(bios[[1]])))]<-NA
  r[rtr_bias(sp1,sp2,bios[[1]])]<-1
  r3<-r
  for(i in 1:iter){
    cat(paste('',i,'-RUN  ',sep=''))
    values(r)<-0
    r[rtr_bias(sp1,sp2,bios[[1]])]<-1
    r2<-r
    r3<-mosaic(r3,r2,fun=sum)
    removeTmpFiles(1) #remove temp files from 1 hour ago - to save memory
  }
  r3[which(is.na(values(bios[[1]])))]<-NA
  r4<-r3
  if(plot==TRUE){
    plot(r3)
  }
  return(r4)}
#### end function