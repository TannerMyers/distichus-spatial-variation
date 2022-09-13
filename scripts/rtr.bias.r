# Custom function from https://github.com/LauraANunes/NunesPearson2017
# From: Nunes, L. A. and Pearson, R. G. (2017), 
# A null biogeographical test for assessing ecological niche evolution.
# Journal of Biogeography. doi: 10.1111/jbi.12910

#### function to store grid cells to test for spatial bias
rtr.bias<-function(sp1,sp2,bios,rotation=TRUE,translation=TRUE){
  library('raster')
  library('maptools')
  library('dismo')
  library('rgeos') #for convexhull and centroid
  sp3<- rbind(unique(sp1),unique(sp2)) # unify both clouds of points to maintain spatial configuration
  poly<-SpatialPoints(sp3) 
  poly<-gConvexHull(poly) #make minimum convex polygon
  cent<-gCentroid(poly) #find centroid of polygon
  
  x_centre<-as.vector(extent(cent))[1]  #find x- centre of polygon
  y_centre<-as.vector(extent(cent))[3]  #find y- centre of polygon
  
  x<-sp3[,1] #x-coordinates
  y<-sp3[,2] #y-coordinates
  
  repeat{
    ###### rotation
    if(rotation==TRUE){ ##option to rotate the points, no translation
      deg<-runif(1,0,360) #random angle in degrees
      
      rad<-deg*0.0174532925 #coverts degrees to radian
      
      R<-matrix(c(cos(rad),sin(rad),0,-sin(rad),cos(rad),0,0,0,1),byrow=T,ncol=3,nrow=3) #clockwise
      
      a<-matrix(c(1,0,0,0,1,0,x_centre,y_centre,1),nrow=3,ncol=3,byrow=T) #rotation at the centre of origin
      
      cn<-matrix(c(1,0,0,0,1,0,-x_centre,-y_centre,1),nrow=3,ncol=3,byrow=T) #rotation at the centre of origin
      
      R<-cn%*%R%*%a
      
      x<-sp3[,1] #x-coordinates
      y<-sp3[,2] #y-coordinates
      
      rot<-matrix(rep(NA),ncol=3,nrow=length(x))
      for(j in 1:length(x)){
        rot[j,]<-R%*%matrix(c(x[j],y[j],1),ncol=1,nrow=3)    #appy rotation to points    
      }
      
      #apply rotation matrix to the original matrix
      repmat<-matrix(rep(c(x_centre,y_centre)),ncol=2,nrow=nrow(sp3),byrow=T)
      
      s<-sp3-repmat
      s<-cbind(s,1)
      s<-matrix(unlist(s),ncol=3,byrow=F)
      so<-tcrossprod(s,round(R,5))
      
      newpool1<-so[,1:2]+repmat     
      
      x<-newpool1[,1]
      y<-newpool1[,2]
    }
    if(translation==TRUE){ #option to translate the points, no rotation
      #### translation 
      xmin=extent(bios)[1];xmax=extent(bios)[2];ymin=extent(bios)[3];ymax=extent(bios)[4]
      
      x_trans<-runif(1,-(xmax-xmin),(xmax-xmin)) #random long trans
      y_trans<-runif(1,-(ymax-ymin),(ymax-ymin)) #random lat trans
      
      x<-round(x+x_trans) #add translation vector
      
      y<-round(y+y_trans) #add translation vector
      
      
      newpool1<-cbind(x,y) #combine new coordinates
    }
    
    # plot(bios[[1]])
    #points(newpool1)
    
    climate_niche = extract(bios, newpool1) #if NA, then a point is in the ocean/outside study region - no variables
    if(isTRUE(length(which(is.na(climate_niche)==FALSE))==length(climate_niche))){  #if TRUE no NAs then break loop, keep RTR polygon
      break}}
  
  gridcells<-extract(bios, newpool1,cellnumbers=TRUE)[,1] #store gridcells to check for duplicates or spatial bias
  
  return(as.vector(gridcells))} #return the gridcells occupied by the simulated points
#### end function ####