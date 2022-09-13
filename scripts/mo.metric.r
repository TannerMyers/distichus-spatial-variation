# Custom function from https://github.com/LauraANunes/NunesPearson2017
# From: Nunes, L. A. and Pearson, R. G. (2017), 
# A null biogeographical test for assessing ecological niche evolution.
# Journal of Biogeography. doi: 10.1111/jbi.12910

########## COMPARE NICHE
#### new metric - niche overlap function ###
mo.metric<-function(sp1,sp2,bios){ 
  
  e_sp1<-na.omit(raster::extract(bios,sp1)) #extract climatic variables for species 1
  e_sp2<-na.omit(raster::extract(bios,sp2)) #extract climatic variables for species 1
  
  e1_max<-apply(e_sp1,2,max) #max of species 1
  e1_min<-apply(e_sp1,2,min) #min of species 1
  
  e2_max<-apply(e_sp2,2,max) #max of species 2
  e2_min<-apply(e_sp2,2,min) #min of species 2
  
  overlap<-rep(NA,ncol(e_sp1)) #to store  overlap of enviromental varialble (Axis overlap, Fig. 1)
  for(n in 1:ncol(e_sp1)){    
    overlap[n]<-(min(e1_max[n],e2_max[n])-max(e1_min[n],e2_min[n]))/(max(e1_max[n],e2_max[n])-min(e1_min[n],e2_min[n]))
    if(overlap[n]<0){overlap[n]<-0} #there is no overlap between the axis if overlap is negative
  }
  
  sumoverlap<-sum(overlap)/ncol(e_sp1) # Cumulative sum of axis overlap, relative to potential overal
  return(list(overlap,sumoverlap)) 
}


################# end function ############