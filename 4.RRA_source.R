# Script for replication of Beyene et al., 2025 
# Author: Menilek Beyene
# Date: April 10, 2025

# Run these in sequence 
source("RRA_IBM.R")
source("RRA_distCC.R")
source("LU_Prop.R")

### Or load data from disk
# RRA_IBM <- read.csv("data/1.raw/2.working/RRA-IBM.csv")
# LU_500<-read.csv("data/1.raw/2.working/RRA-LU500m.csv")

dt <- merge(RRA_IBM, LU_500, by.x = "rra_id", by.y = "id")
rownames(dt) <- dt$rra_id

# Create Land Use Data
dt %>% select(rra_id, x, y, timp,
              Successional.Forest, Meadow, Forest,
              Medium.Density.Residential, High.Density.Residential, Estate.Residential, Rural.Residential,
              Roads, Commercial, Railway, Vacant.Land, Institutional, Industrial,
              Lacustrine, Wetland, Riverine, Beach.Bluff,
              Recreational.Open.Space, Agriculture, Golf.Course) -> RRA_LU

# Create Rapid Restoration Assessment data set
dt %>% select(rra_id, x, y, site_name_RRA, #Site id info
              Project.Manager..1,Surveyor.Name.s.., #Surveyors and Management of RRA
              Date., Assessment.Year., planting_date_RRA, planting_date_IBM, #Date(s) and measure of similarity of planting dates
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
) -> RRA_data

# Create species matrix of planted vegetation in restorations
dt  %>% select(Abies_balsamea:Zizia_aurea) -> RRA_sppMatrix

# # Write to Disk
# write.csv(RRA_data, file = paste0(getwd(), "data/1.raw/2.working/RRA_data.csv")
# write.csv(RRA_LU, file = "RRA_LU.csv")
# write.csv(RRA_sppMatrix, file = "RRA_sppmatrix.csv")