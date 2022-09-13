# Custom function from https://github.com/LauraANunes/NunesPearson2017
# From: Nunes, L. A. and Pearson, R. G. (2017), 
# A null biogeographical test for assessing ecological niche evolution.
# Journal of Biogeography. doi: 10.1111/jbi.12910

################### PLOT SIGNIFICANCE
RTRsignificance<-function(sp1,sp2,bios,rtrs, tails=TRUE, divergence=TRUE){
  
  climate_sp1 = na.omit(extract(bios, sp1)) #extract niches
  climate_sp2 = na.omit(extract(bios, sp2)) #extract niches
  
  observed<-mo.metric(sp1,sp2,bios)[[2]] #observed niche overlap
  meaniche<-mean(rtrs)
 
  ### is PNC, one-tailed test
  
  nineperc<-as.numeric(quantile(rtrs,0.95)) #define critical value limits for up
  top<-observed>as.numeric(quantile(rtrs,0.95)) #niche conserved
  Fn<-ecdf(as.numeric(rtrs))
  location<-Fn(observed)  
  pvalue<-min((sum(rtrs <= observed)/length(rtrs)),(sum(rtrs >= observed)/length(rtrs))) #one tailed each time
  bottom<-NA
  fiveperc<-NA
  
  if(tails==TRUE){ ##two-tailed test
  fiveperc<-as.numeric(quantile(rtrs,0.025)) #define critical value limits for bottom
  nineperc<-as.numeric(quantile(rtrs,0.975)) #define critical value limits for up
  bottom<-observed<as.numeric(quantile(rtrs,0.025)) #niche diverge
  top<-observed>as.numeric(quantile(rtrs,0.975)) #niche conserved
  Fn<-ecdf(as.numeric(rtrs))
  location<-Fn(observed)
  pvalue<-min((sum(rtrs <= observed)/length(rtrs)),(sum(rtrs >= observed)/length(rtrs)))*2 #two tailed p-value
}
  
  if(divergence==TRUE){ ## is PND, one tailed test
    fiveperc<-as.numeric(quantile(rtrs,0.05)) #define critical value limits for bottom
    bottom<-observed<as.numeric(quantile(rtrs,0.05)) #niche diverge
    Fn<-ecdf(as.numeric(rtrs))
    location<-Fn(observed)
    pvalue<-min((sum(rtrs <= observed)/length(rtrs)),(sum(rtrs >= observed)/length(rtrs))) #one tailed
  top<-NA
  nineperc<-NA
  }
   
  hist(rtrs, xlab='Niche Overlap Value', ylab='Frequency',main=NULL,xlim=c(0,1))
  legend('topright',legend=c('Observed Niche Overlap Value'), lty=c(1),lwd=c(2),col='red')
  abline(v=observed, col='red', lwd=2)
  return(data.frame(observed,bottom,top,meaniche,fiveperc,nineperc,pvalue,location))}

################## end function