# Script for replication of Beyene et al., 2025 
# Author: Menilek Beyene
# Date: April 10, 2025

# Fuzzy Text Matching in R -----------------------------------------------------

# Objective 
"The objective of this script to is to combine two datasets via a fuzzy matching. 
Fuzzy matching is completed via the fedmatch functions. 
This involves a data cleaning step followed by a similarity determination."

# Library and Packages ---------------------------------------------------------
library(fedmatch)
library(data.table)
library(tidyverse)

# Step 1 - Clean data set of errors or typos -----------------------------------
# Load data
RRA<-read.csv("data/1.raw/RRA_data.csv")
IBM<-read.csv("data/1.raw/IBM_data.csv")
# Rename columns for site names in RRA and IBM to assist matching 
RRA %>% dplyr::rename(site_name = Site.Name.) -> RRA 
IBM %>% dplyr::rename(site_name = Site.Name) -> IBM 
# To replace special cases "_" with spaces " " 
IBM$site_name <- stringr::str_replace_all(IBM$site_name, "_", " ")
  # OR use FedMatch::clean_strings() default
  # https://cran.r-project.org/web/packages/fedmatch/vignettes/Intro-to-fedmatch.html 

# Replace site name entries with cleaned text
RRA$site_name <- clean_strings(RRA$site_name)
IBM$site_name <- clean_strings(IBM$site_name)

# Additional cleaning of sites names by takeing string input (i.e., column of names)
clean_sitenames <- function(sitenames) {
  # sitenames <- tolower(sitenames) #convert all to lowercase
  # sitenames <- stringr::str_replace_all(sitenames,"[^a-zA-Z0-9\\s]", "") #replace all non-text or numbers from string
  sitenames <- stringr::str_remove_all(sitenames, "compensation")
  sitenames <- stringr::str_remove_all(sitenames, "comp")#remove instances of compensation
  sitenames <- trimws(sitenames, "both") #remove leading and trailing whitespace

  return(sitenames) #output cleaned site names
}

# Remove instances of compensation or comp from words
IBM$site_name <- clean_sitenames(IBM$site_name)
RRA$site_name <- clean_sitenames(RRA$site_name)

# Remove rows that don't meet the count check criteria (i.e., NA)
IBM <- IBM[IBM$count_check != "#N/A",] #remove rows with NA in the check column

# Fill rows with NA is species matrix with zeros
IBM %>% replace(is.na(.), 0) -> IBM

# produce a unique key for each entry
RRA$rra_id <- seq(1:nrow(RRA))

# To help only work with columns of interest
RRA %>% dplyr::select(rra_id, site_name, x, y, #Site id info
               Survey_team, #Surveyors and Management of RRA
               Date.,Assessment.Year.,year, #Date(s)
               Applicable.project.deliverable.categories., Applicable.Project.Deliverable.Categories., #Ecosystem deliverable
               Natural.cover.component.score., Which.planting.types.are.a.part.of.the.deliverable., #Scoring of total natural cover
               Survival.score., Average.health.of.surviving.plants., ConifVegScore, #Coniferous plantings
               Survival.score..1, Average.health.of.surviving.plants..1, DecidVegScore, #Deciduous plantings
               Survival.score..2, Average.health.of.surviving.plants..2, ShrubVegScore, #Shrub plantings
               Survival.score..3, Average.health.of.surviving.plants..3, HerbVegScore, #Herbaceous plantings
               Survival.score..4, Average.health.of.surviving.plants..4, BioEngVegScore, #BioEngineering plantings
               Survival.score..5, Average.health.of.surviving.plants..5, CaliperVegScore, #Caliper plantings - May be in reference to caliper trees (or larger single trees i.e., Ufor)
               Invasive.species., Extent.of.invasive.species., Invasive.Score., Aggressive.invasive.undesireable.herbaceous.coverage....., #Invasive species data
               Year.last.treated., Invasive.species.observed.in.treated.location., #Invasive species management
               Native.species., Natural.native.regeneration.score., Is.native.regeneration.of.a.type.consistent.with.the.deliverable.goals.present., Native.wildfolwer.and.grass.coverage....., #Native species natural regeneration
               Browsing.of.planted.vegitation., Extent.of.browsing., Browse.Score., Select.sources.of.browse.that.have.an.impact.on.project.success., Type.of.insects., #Browsing magnitude
               Installations.present.to.deter.browse., Select.the.browse.deterrents.needing.repair..replacement..installation., Notes.on.browse.deterrents.needed., #Browsing management
               Observed.challenges.to.plant.survival., Infill.recommendations., Notes.on.infill.recommendations., #Notes and interpretation
) -> RRA_sites
# To help work with only columns of interest
IBM %>% dplyr::select(IBM_id, site_name, day_month, year, 
               Entry_by:last_col()
) -> IBM_sites

# Need to create a planting date column for the RRA
# Convert all dates to ISO standards (yyyy-mm-dd)
RRA_sites$Date. <- str_replace_all(RRA_sites$Date., ", "," ")
RRA_sites$Date. <- str_replace_all(RRA_sites$Date., " ", "/")

#Convert to a date format for following step
RRA_sites$Date. <- as.Date(RRA_sites$Date., format = "%b/%d/%Y") #convert to date format

# IBM_sites$day_month <- as.Date(IBM_sites$day_month, format = "%d-%b")
IBM_sites %>% dplyr::rename(planting_year = year) -> IBM_sites
IBM_sites$planting_date <- paste(IBM_sites$day_month, IBM_sites$planting_year, sep = "-")  
IBM_sites$planting_date<- as.numeric(as.Date(IBM_sites$planting_date, format = "%d-%b-%Y"))

# Assessment date minus the number of years since planting should equal the planting date
RRA_sites$planting_date <- as.numeric(as.Date(ymd(RRA_sites$Date.) - years(RRA_sites$Assessment.Year.))) #convert date to numeric to allow rescaling during difference distance measure





# Step 2 - Fuzzy Matches -------------------------------------------------------
# Match site_names in RRA and IBM data sets with RRA as the primary key

##Seminal Paper on record linkage
#https://www.tandfonline.com/doi/pdf/10.1080/01621459.1969.10501049

# Complete multivar matches for each observation in RRA to determine closest match in data2
result1 <- fedmatch::merge_plus(
  data1 = RRA_sites,
  data2 = IBM_sites, 
  #Multi variable record linkage
  match_type = "multivar",
  by = c("site_name", "planting_date"),
  suffixes = c("_RRA", "_IBM"), 
  unique_key_1 = "rra_id", unique_key_2 = "IBM_id",
  allow.cartesian = T,
  multivar_settings = build_multivar_settings(
    compare_type = c("wgt_jaccard_dist", "difference"),
    wgts = c(.8, .1)
  )
)

RRA_IBM <- as.data.frame(result1$matches)
RRA_IBM %>% dplyr::select(rra_id, IBM_id, x, y, 
                         site_name_RRA, site_name_IBM, 
                         Date., Assessment.Year., planting_date_RRA, 
                         planting_date_IBM, day_month,
                         Applicable.Project.Deliverable.Categories., 
                         Applicable.project.deliverable.categories., 
                         Which.planting.types.are.a.part.of.the.deliverable.) -> compare_matches #for viewing matches


compare_matches$planting_date_RRA <- as.Date(compare_matches$planting_date_RRA)
compare_matches$planting_date_IBM <- as.Date(compare_matches$planting_date_IBM)


# Fuzzy Matching ---------------------------------------------------------------
# # Complete fuzzy matches for each observation in RRA to determine closest match in data2
# result1 <- fedmatch::merge_plus(
#   data1 = RRA_sites,
#   data2 = IBM_sites, 
#   
#   #Multi variable record linkage
#   match_type = "fuzzy",
#   by = c("site_name"),
#   suffixes = c("_RRA", "_IBM"), 
#   unique_key_1 = "rra_id", unique_key_2 = "IBM_id",
#   allow.cartesian = T,
#   fuzzy_settings = build_fuzzy_settings(method = "wgt_jaccard", nthread = 2, 
#                                         maxDist = 0.5,
#     # wgts = c(0.8),
#     # maxDist = 
#     )
#   )
# 
# RRA_IBM <- as.data.frame(result1$matches)
# RRA_IBM %>% select(IBM_id, rra_id, x, y, 
#                          site_name_IBM, site_name_RRA,
#                          # Date., Assessment.Year., planting_date_RRA, planting_date_IBM, day_month,
#                          # Applicable.Project.Deliverable.Categories., Applicable.project.deliverable.categories., Which.planting.types.are.a.part.of.the.deliverable.
#                          ) -> compare_matches #for viewing matches

# Tier Matching ----------------------------------------------------------------
# tier_list <- list(
#   a = build_tier(match_type = "exact"),
#   b = build_tier(match_type = "fuzzy", 
#                  fuzzy_settings = build_fuzzy_settings(method = "wgt_jaccard", nthread = 2, 
#                                                        maxDist = 0.8)),
#   c = build_tier(match_type = "multivar", multivar_settings = build_multivar_settings(
#     compare_type = c("wgt_jaccard_dist", "difference"),
#     wgts = c(.8, .1), 
#     nthread = 2
#     # logit = NULL, missing = FALSE, wgts = 1,
#     # compare_type = "stringdist", blocks = NULL, blocks.x = NULL, blocks.y = NULL,
#     # top = 1, threshold = NULL
#   ))
#   
# )
# 
# 
# 
# 
# tier_results <- tier_match(data1 = IBM_sites,
#                            data2 = RRA_sites,
#                            by = c("site_name"),
#                            
#                            suffixes = c("_IBM", "_RRA"), 
#                            unique_key_1 = "IBM_id", unique_key_2 = "rra_id",
#                            tiers = tier_list, takeout = "neither", verbose = TRUE,
#                            score_settings = build_score_settings(score_var_x = "IBM_sites",
#                                                                  score_var_y = "RRA_sites",
#                                                                  wgts = 0.8,
#                                                                  score_type = c("stringdist", "difference")
#                                                                  )
#                            )



# ----------------------------------------------------------
# Why do we have the timp merge down below? Investigate this.  

### Write files to working directory ###
# If you want to create heard copies this can be done using these lines
# These write a full list of files but we want 
# write.csv(RRA_IBM, file = "data/1.raw/2.working/RRA_IBM2.csv") #Full matched records
# write.csv(compare_matches, file = "data/1.raw/2.working/comparison_matches.csv") #Comparison records
# write.csv(IBM_sitenames, file = "IBM_sites.csv")