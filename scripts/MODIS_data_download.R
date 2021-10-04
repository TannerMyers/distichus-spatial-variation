setwd("~/Dropbox/Distichus_Project/distichus_spatial_morphological_variation/")

library(MODISTools)

products <- mt_products()
products
  # MOD13Q1 is the MODIS product we need

mt_bands("MOD13Q1")
# band                        description
# 1           250m_16_days_blue_reflectance         Surface Reflectance Band 3
# 2  250m_16_days_composite_day_of_the_year               Day of year VI pixel
# 3                        250m_16_days_EVI                 16 day EVI average
# 4            250m_16_days_MIR_reflectance         Surface Reflectance Band 7
# 5                       250m_16_days_NDVI                16 day NDVI average
# 6            250m_16_days_NIR_reflectance         Surface Reflectance Band 2
# 7          250m_16_days_pixel_reliability    Quality reliability of VI pixel
# 8            250m_16_days_red_reflectance         Surface Reflectance Band 1
# 9     250m_16_days_relative_azimuth_angle Relative azimuth angle of VI pixel
# 10          250m_16_days_sun_zenith_angle       Sun zenith angle of VI pixel
# 11         250m_16_days_view_zenith_angle      View zenith angle of VI Pixel
# 12                250m_16_days_VI_Quality              VI quality indicators

# Load dataframe with sampling locality coordinates
data <- read.csv("data-coords.csv",header = TRUE)
data <- data[,1:2]

dates <- mt_dates(product = "MOD13Q1", lat = 19, lon = -72)
head(dates)

  
mt_batch_subset(df= data, product = "MOD13Q1",
                band ="250m_16_days_EVI",
                start = "2000-02-18",
                end = "2020-10-31",
                out_dir="environmental-data/MODIS")

mt_batch_subset(df= data, product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                start = "2000-02-18",
                end = "2020-10-31",
                out_dir="environmental-data/MODIS")

