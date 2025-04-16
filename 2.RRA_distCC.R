#### Distance to Single City Center ####


### Read set set of points and determine there distance to a central single point
city_filepath <- paste0(getwd(), "/data/1.raw/city_halls_latlong.csv")
point_filepath <- paste0(getwd(),"/data/1.raw/RRA_data.csv")

city_halls <-  read.csv(city_filepath)
points <- read.csv(point_filepath)
epsg <- as.numeric("4326")

city_halls %>% filter(Name == "Toronto City Hall") -> TO_city_hall


# load point csv
point_data <- read.csv(point_filepath)
# Create unique RRA_id numbers
point_data$id <- seq(nrow(point_data))


point_data %>% dplyr::select(id, x, y) -> points_sel

##convert csv to sf class files
# convert to simple feature POINT 
RRA_distCC <- st_as_sf(points_sel, coords = c("x", "y"), crs = epsg)
TO_city_hall <- st_as_sf(TO_city_hall, coords = c("x", "y"), crs = epsg)

#Change projection

###Calculate distances
points$dist_TO_cityhall <- st_distance(point_sel, TO_city_hall)
points_sel$dist_TO_cityhall <- st_distance(point_sel, TO_city_hall)

rm(city_halls, points, point_sel, TO_city_hall)
# # Write distance in meters from city center for each RRA to disk
# write.csv(points_sel, "/data/1.raw/2.working/RRA_distCC.csv")