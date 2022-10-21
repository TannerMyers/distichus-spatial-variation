
library(ENMTools)

# Load filtered occurrences downloaded with `spocc` R package
occs_ad <- read_csv(paste0(occ_dir, "filtered_distichus_occurrences.csv"))

# Add imputed rasters to a new raster stack
env_imputed <- raster::stack(list.files(path = "data/env_imputed", pattern = ".asc", full.names = TRUE))

# Subset rasters left after tests for collinearity
env2 <- env_imputed[[c("CHELSA_bio10_03", "CHELSA_bio10_04", "CHELSA_bio10_05",
            "CHELSA_bio10_15", "CHELSA_bio10_16",
            "march_EVI_mean", "may_NDVI_mean")]]