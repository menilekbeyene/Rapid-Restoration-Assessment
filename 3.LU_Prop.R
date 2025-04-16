# Script for replication of Beyene et al., 2025 
# Author: Nicholas Sookhan & Menilek Beyene
# Date: April 10, 2025


#######################################################################################################
### GET COMMAND LINE ARGUMENTS ### 
#######################################################################################################

# if command line arguments exist, then get args from command line
if(!is.na(commandArgs(trailingOnly = T)[1])) { 
  # full path to site list
  point_path <- commandArgs(trailingOnly = T)[1]
  # full path to land cover raster
  lc_path <- commandArgs(trailingOnly = T)[2]
  # projection of site list csv
  epsg <- as.numeric(commandArgs(trailingOnly = T)[3])
  # buffer size
  bsize <- as.numeric(commandArgs(trailingOnly = T)[4])
  # X column
  X <- commandArgs(trailingOnly = T)[5]
  # Y column
  Y <- commandArgs(trailingOnly = T)[6]
  # ID 
  id <- commandArgs(trailingOnly = T)[7]
  
  # else if arguments already exist then do not redefine
} else if (exists("point_path") & exists("lc_path") & exists("epsg") & exists("bsize") & exists("X") & exists("Y") & exists("id")) {
  
  # define arguments here 
} else {
  # full path to site list
  # point_path <- paste0(getwd(),"/data/1.raw/RRA_data.csv")
  # full path to land use raster (set to where raster has been saved)
  # GeoJson URL
  lc_path <- paste0(getwd(), "/data/1.raw/TRCA_LU-2017_epsg4326.tif")
  # projection of site list csv
  epsg <- as.numeric("4326")
  # buffer size
  bsize <- as.numeric("500") #adjust buffer size as needed to measure at different scales (25, 50, 100, 200)
  # X
  X <- "x"
  # Y 
  Y <- "y"
  # id 
  id <- "id"
}

#######################################################################################################
### LOAD PACKAGES ###
#######################################################################################################
ip <- installed.packages()[,1]
libs <- c("sf","raster","parallel","doParallel","data.table","stringr","rgdal")
# install packages that are not installed
if (sum(!(libs %in% ip))>0){ install.packages(libs[!(libs %in% ip)],repos="https://cloud.r-project.org/") }
# load packages
lapply(libs,library, character.only=T)
library(landscapemetrics)

#######################################################################################################
### GLOBAL OPTIONS ###
#######################################################################################################
rasterOptions(maxmemory = 1e+10)
options(stringsAsFactors = F)

#######################################################################################################
### LOAD DATA ###
#######################################################################################################
source("UFRE_distCC.R")

# load land cover csv
lc <- raster(lc_path)

#######################################################################################################
### MANAGE DATA ###
#######################################################################################################

### LAND COVER RASTER ### 
# get landcover raster  projection
lc_proj <- proj4string(lc)

### POINTS ###
# only keep relevant columns 
point <- point_data[,c("id","x","y")]
# drop points with missing x and/or missing y
point <- point[!(is.na(point[,"x"])|is.na(point[,"y"])),]
# convert to simple feature POINT 
point <- st_as_sf(point, coords = c("x", "y"), crs = epsg)
# transform to projection of raster
point <- st_transform(point,lc_proj)
# create buffer 
buffer <- st_cast(st_buffer(point,bsize),"POLYGON")
# convert buffers to list 
buffer <- split(buffer[,1], f = buffer[,1,drop=T])

#######################################################################################################
### CALC METRICS ###
#######################################################################################################

### calculate metrics ### 
# use cores - 2 
UseCores <- detectCores() - 2

#Register CoreCluster                                                                                                                                                                                                                                                                         
cl <- makeCluster(UseCores)                 
registerDoParallel(cl)

system.time(
  point_metric <- 
    foreach(i=1:length(buffer)) %dopar% {
      library(raster)
      library(sf)
      library(data.table)
      
      # check if buffer[[i]] overlaps lc
      if ( is.null( intersect(extent(buffer[[i]]), extent(lc))) ) {
        # if NO overlap 
        metric_calc <- NA
        prop_missing <- data.frame("prop_missing" = 1)
      } else {
        # if YES overlap 
        #crops to extent
        clip1 <- crop(lc, extent(buffer[[i]]))                                                              
        #crops to polygon edge & converts to raster
        clip2 <- rasterize(buffer[[i]], clip1, mask=TRUE)
        # if overlap then check if populated cells occur
        if ( cellStats(clip2>0, stat = "sum") == 0 ) {
          # if NO populated cells
          metric_calc <- NA
          prop_missing <- data.frame("prop_missing" = 1)
        } else {
          # create mask 
          clip3 <- rasterize(buffer[[i]], setValues(clip1,1), mask=TRUE)
          # calculate PLAND
          pland = table((matrix(clip2)))/cellStats(clip3,sum)
          # total area by class
          area = table((matrix(clip2))) * (res(lc)[1]*res(lc)[2])
          # total area missing 
          missing = 1 - (cellStats(!is.na(clip2),sum) / cellStats(!is.na(clip3),sum))
          # coerce to dataframe
          pland = matrix(pland,ncol=dim(pland),dimnames = list("1",paste("pland", names(pland),sep="_")))
          area = matrix(area,ncol=dim(area),dimnames = list("1",paste("area", names(area),sep="_")))
          # 
          metric_calc = data.frame(pland, area)
          prop_missing = data.frame("prop_missing" = missing)
        }
      }
      return(list("metric_calc"=metric_calc,"prop_missing"=prop_missing))
    })
#end cluster
stopCluster(cl)

#######################################################################################################
### ORGANIZE METRICS ###
#######################################################################################################
# assign site names to metric list
names(point_metric) <- names(buffer)

### SUMMARY ### 
point_metric_report <-  rbindlist( lapply(point_metric,function(x)x$prop_missing), idcol = T,fill=T)
# change .id to site_id
colnames(point_metric_report)[1] <- id
# get missing
missing <- point_metric_report$site_id[round(point_metric_report$prop_missing*100,2)==100]

### METRIC VALUES ###
point_metric_calc <- lapply(point_metric,function(x)x$metric_calc)
# remove missing 
point_metric_calc <- point_metric_calc[!names(point_metric_calc) %in% missing]
#Some list items are not data.frame data.table or list objects must be removed
  #Issue is due to centroids extending beyond extent of the Land Use raster

#Find na values in list of list of data.tables
which_na<-which(rapply(point_metric_calc, is.na))
#iterative removal of NA entries
for(i in 1:length(which_na)) {
  point_metric_calc[[names(which_na[i])]]=NULL #remove each element that is na from list of list
}
# to data.table
point_metric_calc <- rbindlist(point_metric_calc, idcol = T,fill=T)
setDF(point_metric_calc)
point_metric_calc[is.na(point_metric_calc)] = 0
colnames(point_metric_calc)[1] <- id

# Load conversion table for LU column names
key_path <- paste0(getwd(),"/data/1.raw/TIMP_LU_conversion_table.csv")
# Set land use categories
TIMP_LU <- read.csv(file = key_path)
# Create key by which we will rename colnames in data
TIMP_LU %>% dplyr::select(LU_id | Landuse) -> idkey

dt <- point_metric_calc
# Selecting columns that contain land use propertions and id column
dt %>% dplyr::select( "id" | starts_with("pland")) -> dt_pland 

# #Select columns with area of landuses 
# dt %>% dplyr::select( "id" | starts_with("area_")) -> dt_area

# match dataframe column names to key by removing pland prefix from column names
dt_pland %>% rename_all(~stringr::str_replace(.,"^pland_","")) -> dt_pland  
# dt_pland %>% rename_with( ~ paste0("pland_, .x")) alternate code to rename by adding prefix
#Rename columns in working dataframe
colnames(dt_pland) <- dplyr::recode(
  colnames(dt_pland),
  !!!setNames(as.character(idkey$Landuse), idkey$LU_id)
)


#######################################################################################################
### WRITE TO DISK ### 
#######################################################################################################
output_name_report <- paste0(str_extract(point_path,".+(?=\\.csv)"),"_metric","_",bsize,"_report",".csv")
output_name_calc <- paste0(str_extract(point_path,".+(?=\\.csv)"),"_metric","_",bsize,"_calc",".csv")
output_name_prop <- paste0(str_extract(point_path,".+(?=\\.csv)"),"_metric","_",bsize,"_prop",".csv")
# 
write.csv(point_metric_report, file=output_name_report, row.names = F)
write.csv(point_metric_calc, file=output_name_calc, row.names = F)
write.csv(dt_pland, file=output_name_prop, row.names = F)