#### begin ####
library(MODIS)
library(terra)
library(sf)
library(raster)
# step zero: set year
year <- "2023"

# step 1 : read each year's landcover data
hdf_files <- list.files(path = "data/modis/MCD12Q1_061-20241120_093859/", pattern = paste(year,"001.*\\.hdf$", sep= ""), full.names = TRUE)
output_file <- paste("tifs/MCD12Q1_061-",year,".tif", sep = "")

# step2: read the pair of hfd files

sds_info1 <- sds(hdf_files[3])
hdf2 <- rast(hdf_files[4])
# List subdatasets (individual layers)
subdatasets <- sds(hdf)

sds_info2 <- sds(hdf_files[4])

#step3 : set crs and merge
WG <- st_read("../shapefiles/S_WG.shp")
#layer 1
layer_test <- sds_info1$MCD12Q1.A2023001.h25v07.061.2024252132124
layer_rast <- raster(layer_test)
WG_crs <- st_crs(WG)$wkt
layer_test_proj <- projectRaster(layer_rast, crs = crs(WG))


layer_test2 <- project(layer_test, WG_crs)
crs(layer_test) <- "+proj=longlat +datum=WGS84 +no_defs"          
#writeRaster(layer_test, "test-rastercrs.tif")
#layer 2
layer_best <- sds_info2$MCD12Q1.A2023001.h25v08.061.2024252132154
crs(layer_best) <- "+proj=longlat +datum=WGS84 +no_defs"          
#writeRaster(layer_best, "best-rastercrs.tif")
# mosaic
test_mosaic <- terra::mosaic(hdf$LC_Type2, hdf$LC_Type2)
crs(test_mosaic)
plot(st_geometry(test_mosaic))


writeRaster(test_mosaic, output_file)

#### turn into a for loop ####

for (year in seq(2014,2023)) {
    hdf_files <- list.files(path = "data/modis/MCD12Q1_061-20241120_093859/", pattern = paste(year,"001",".*\\.hdf$", sep= ""), full.names = TRUE)
    print(hdf_files)
    output_file <- paste("tifs/MCD12Q1_061-",as.character(year),".tif", sep = "")
    
    # step2: read the pair of hfd files
    
    hdf1 <- rast(hdf_files[1])
    hdf2 <- rast(hdf_files[2])
    
    layer_test <- hdf1$LC_Type2
    layer_best <- hdf2$LC_Type2

    #writeRaster(layer_test, "test-rastercrs.tif")
    #layer 2
    #writeRaster(layer_best, "best-rastercrs.tif")
    
    # step 3: mosaic
    test_mosaic <- terra::mosaic(layer_test, layer_best)
    test_mosaic
    writeRaster(test_mosaic, output_file, overwrite = TRUE)
    
}
