## Load function for imputing NAs in the raster from
## https://stackoverflow.com/questions/45641168/fill-in-gaps-e-g-not-single-cells-of-na-values-in-raster-using-a-neighborhood
fill.na <- function(x) {
  center = 0.5 + (width*width/2)
  if(is.na(x)[center]) {
    return(mean(x, na.rm = TRUE))
  } else {
    return(x[center])
  }
}