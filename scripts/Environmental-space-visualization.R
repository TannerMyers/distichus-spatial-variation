setwd("~/Dropbox/Distichus_Project/distichus_spatial_morphological_variation/")

# Load data for individuals that have been measured 
data <- read.csv("complete-data.csv", header=TRUE)
columns <- c("Subspecies", "Longitude", "Latitude", "Color")

df <- data[, columns]

# Load processed MERRAClim bioclimatic variables at 2.5 arc seconds resolution
bioclim <- stack(list.files(path="merra_new", pattern = ".asc", full.names = T))

# Extract raster values for the localities in `df`
bc_values <- raster::extract(bioclim,data.frame(df$Longitude, df$Latitude))

# Combine into one dataframe
df <- cbind(df, bc_values)

leg.txt <- c(expression(italic("ignigularis")),
             expression(italic("ravitergum")),
             expression(italic("favillarum")),
             expression(italic("properus")),
             expression(italic("dominicensis")),
             expression(italic("aurifer")),
             expression(italic("suppar")),
             expression(italic("vinosus")))

pdf("MERRAclim-environmental-biplot.pdf")
plot(df$bio12 ~ df$bio1, 
     col=df$Color, pch=16, cex=1.25,
     xlab="Annual Mean Temperature", ylab="Annual Precipitation",
     main="Environmental Space (MERRAClim)")
legend("topleft", leg.txt, col=unique(data$Color), pch =16,
       cex=1)
dev.off()


data <- read.csv("complete-data.csv", header=TRUE)
  head(data)
  str(data)

  unique(data$Subspecies)
# Levels: aurifer dominicensis favillarum ignigularis properus ravitergum suppar vinosus
  
  data$Color[data$Subspecies=="aurifer"] <-"peru" 
  data$Color[data$Subspecies=="vinosus"] <-"firebrick3" 
  data$Color[data$Subspecies=="suppar"] <-"purple1" 
  data$Color[data$Subspecies=="favillarum"] <-"dodgerblue" 
  data$Color[data$Subspecies=="ravitergum"] <-"yellow" 
  data$Color[data$Subspecies=="dominicensis"] <-"bisque" 
  data$Color[data$Subspecies=="ignigularis"] <-"darkorange2" 
  data$Color[data$Subspecies=="properus"] <-"seagreen3"   

  leg.txt <- c(expression(italic("ignigularis")),
               expression(italic("ravitergum")),
               expression(italic("favillarum")),
               expression(italic("properus")),
               expression(italic("dominicensis")),
               expression(italic("aurifer")),
               expression(italic("suppar")),
               expression(italic("vinosus")))
  
  pdf("Worldclim-environmental-biplot.pdf")
  plot(data$bio12_23 ~ data$bio1_23, 
       col=data$Color, pch=16, cex=1.25,
       xlab="Annual Mean Temperature", ylab="Annual Precipitation",
       main="Environmental Space (WorldClim)")
  legend("topleft", leg.txt, col=unique(data$Color), pch =16,
         cex=1)
  dev.off()
  