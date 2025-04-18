# Analysis of Rapid Restoration Assessment Scores ####
# Script for replication of Beyene et al., 2025 
# Author: Menilek Beyene
# Date: April 10, 2025




#### Load Library of Packages ####
library(tidyverse)
library(fedmatch)
library(data.table)
library(vegan)
library(ape)
library(dplyr)
library(ggcorrplot)
library(ordinal) #What does clmm use? Maximum likelihood vs restricted maximum likelihood
library(sjlabelled) #Set labels for columns
library(car)

# library(emmeans)
# # library(RVAideMemoire)
################# Figs, Tables, and Supplementary Materials ###################
# library(ggcorrplot)
# library(olsrr)
# library(MASS)
# library(stargazer)

# # If loading from disk import using lines below
# #### Import data ####
RRA_IBM <- read.csv("data/1.raw/2.working/RRA_data.csv") #Rapid Restoration Assessment (RRA) and Internal Billing Memo's (IBM) data
LU_path <- "data/1.raw/2.working/RRA_LU.csv"
RRA_LU <- read.csv(LU_path) #Land use percentages at a 500 m radii
RRA_distCC <- read.csv("data/1.raw/2.working/RRA_distCC.csv")


# # If calling source
# source("4.RRA_source.R")
# 
# head(RRA_LU)
# head(RRA_data)
# head(RRA_sppMatrix)
# head(RRA_distCC)

RRA_IBM <- RRA_IBM %>% rename(RRA_id = "rra_id")
RRA_distCC <- RRA_distCC  %>% rename(RRA_id = "id") 
#### Transform Data ####
RRA_LU %>% filter(rra_id %in% RRA_IBM$RRA_id) %>%
  dplyr::rename(RRA_id = rra_id) -> RRA_LU_filt

### Join data sets ###
RRA_data <- RRA_IBM %>% left_join(RRA_LU_filt, by = "RRA_id") #Join land use dataset to the combined RRA-Species matrix data

#Add the distance to city center (Toronto city hall) as a column in combined data
RRA_data$dist_TO <- RRA_distCC %>%
  dplyr::filter(RRA_id %in% 
                  RRA_data$RRA_id) %>% dplyr::group_by(RRA_id) %>% 
  dplyr::summarise(dist_TO_cityhall = mean(dist_TO_cityhall)) %>% 
  dplyr::select(dist_TO_cityhall) %>% unlist() %>% as.numeric() #Create distance to downtown city hall - urbanization metric

#### Check data classes and structure ####
### Set levels to increase towards positive outcome (i.e., 0->1->2->3->4) for more intuitive understanding
'Before fitting an ordinal logistic regression, check the ordering of the outcome 
variable. An OR > 1 corresponds to a risk factor that is associated with greater 
probability of higher levels of the outcome variable. Therefore, typically, the 
ordering desired for an outcome variable is from less severe to more severe. 
Whatever the ordering desired, check using levels() and, if needed, reorder using factor().'

RRA_data$site_name <- as.factor(RRA_data$site_name_RRA) #Set site names as factors
#Convert response variables of restoration outcomes to ordinal - survival, health, regeneration, etc. 
RRA_data[,c(37, 39:56, 67)] <- lapply(RRA_data[,c(37, 39:56, 67)], function(x) {
  #Convert columns to ordered factors
  ordered(x, sort(unique(x))) #To set desired direction the original order is reversed for SS, HASP,
}
)

#Convert response variables of disturbance outcomes to ordinal - herbivory, humans, invasive species, etc. 
RRA_data[,c(92:97, 103:105)] <- lapply(RRA_data[,c(92:97, 103:105)], function(x) {
  #Convert columns to ordered factors
  ordered(x, rev(sort(unique(x)))) #reverse order of factors to increase towards positive outcome!***
}
)


#### Species Composition Analysis using PCoA ####
##Hellenger transformation - https://r.qcbs.ca/workshop09/book-en/transformations.html#hellinger-transformation 
rownames(RRA_data) <- RRA_data$RRA_id #Set row names to site IDs 
# rownames(spe) <- RRA_data$RRA_id
spe <- RRA_data %>% filter(Applicable.project.deliverable.categories. != "Meadow" &
                             Applicable.project.deliverable.categories. != "Green Infrastructure" &
                             Applicable.project.deliverable.categories. != "Shoreline") %>% 
  dplyr::select(RRA_id, Applicable.project.deliverable.categories., 
                Abies_balsamea:Zizia_aurea, -Salix_species) #Create species matrix and remove hybrid salix species (heavily weighs one axis)

spe <- as_tibble(spe %>% dplyr::select(Abies_balsamea:Zizia_aurea)) #Create species matrix and remove hybrid salix species (heavily weighs one axis)
spe.code <- str_to_title(paste0(str_split_i(as.character(colnames(spe)), "_", i = 1), " ",
                                str_split_i(as.character(colnames(spe)), "_", i = 2)))  #The brackets around the period avoid the spatial character nature of period in regular expressions which means "any single character"
colnames(spe) <- gsub("(\\b[A-Z][a-z][a-z])|.", "\\1", spe.code, perl=TRUE) #Change column names to species codes
spe <- tibble(spe, .name_repair = "unique") #Assign unique names to each column to avoid combining species

spe <- as.data.frame(spe)


spe.dt <- sapply(spe, as.numeric) %>% as.data.frame()
spe.dt$RRA_id <- RRA_data %>% filter(Applicable.project.deliverable.categories. != "Meadow" &
                                          Applicable.project.deliverable.categories. != "Green Infrastructure" &
                                          Applicable.project.deliverable.categories. != "Shoreline") %>% dplyr::select(RRA_id) %>% as.vector()

spe.dt <- spe.dt %>% filter(rowSums(across(AbiBal:ZizAur)) !=0) %>% #Remove empty rows
  # column_to_rownames('RRA_id') %>% #Convert id column to rownames to retain rownames in filtered data
  dplyr::select(AbiBal:ZizAur) #Select species columns

colnames(spe.dt)
rownames(spe.dt) 
nrow(spe) > nrow(spe.dt) #Check if empty rows were removed (should be lower then )

spe.dt <- spe.dt[,colSums(spe.dt) !=0] # check shows that spe has no colSums == 0 

library(ade4)
library(vegan)
library(ape)
library(gclus)
#Species abundances
library(BiodiversityR)
spe.abund <- rowSums(as.data.frame(do.call(rbind,lapply(spe[,-c(1,2,258)],as.numeric))))

#Here I work through a question about whether to pre-transform abundance data...
#The basis is that pre-transformation would result in euclidean distances. 
spe.dist <- vegdist(spe.dt, method = "bray") #bray-curtis dissimilarity used for quantitative species abundances
spe.hel <- decostand(spe.dt[,c(1:126)], method = "hellinger") ##Control for values that are very large (a problem for some sites)
spe.dist2 <- vegdist(spe.hel) #Eucliean dissimilarity used for quantitative species abundances
spe.dist3 <- vegdist(spe.dt) #Euclidean dissimilarity used for quantitative species abundances

# Ordination of site x species abundance matrix 
spe.cmd <- cmdscale(spe.dist, #cmdscale (Classical Multidimensional Scaling) produces a PCoA or spread of points equal to their (dis)similarity. 
                      k=3, # Maximum dimensions of space the data should be represented in 
                      # Three was selected to capture functional group differences and reduce number of new parameters in later regression. 
                      eig = T, # We want to return eigenvalues. 
                      add = T) # Add = T computes a minimal additive constant c* that makes dissimilarities Euclidean* Important for later clustering. 

#This should be correct as species abundances have been transformed to account for abundance differences.
spe.h.cmd <- cmdscale(spe.dist2, #cmdscale (Classical Multidimensional Scaling) produces a PCoA or spread of points equal to their (dis)similarity.
                    k=3, # Maximum dimensions of space the data should be represented in
                    # Three was selected to capture functional group differences and reduce number of new parameters in later regression.
                    eig = T, # We want to return eigenvalues.
                    add = T) # Add = T computes a minimal additive constant c* that makes dissimilarities Euclidean* Important for later clustering.

#This should be incorrect as species abundances have not been transformed to account for abundance differences.
spe.cmd.ne <- cmdscale(spe.dist3, #cmdscale (Classical Multidimensional Scaling) produces a PCoA or spread of points equal to their (dis)similarity.
                    k=3, # Maximum dimensions of space the data should be represented in
                    # Three was selected to capture functional group differences and reduce number of new parameters in later regression.
                    eig = T, # We want to return eigenvalues.
                    add = T) # Add = T computes a minimal additive constant c* that makes dissimilarities Euclidean* Important for later clustering.
# spe.h.cmd2 <- cmdscale(spe.dist4, #cmdscale (Classical Multidimensional Scaling) produces a PCoA or spread of points equal to their (dis)similarity. 
#                        k=3, # Maximum dimensions of space the data should be represented in 
#                        # Three was selected to capture functional group differences and reduce number of new parameters in later regression. 
#                        eig = T, # We want to return eigenvalues. 
#                        add = T) # Add = T computes a minimal additive constant c* that makes dissimilarities Euclidean* Important for later clustering. 

### WAScores

spe.cmd$species <- wascores(spe.cmd$points, # Ordination scores for sites.
                              spe.dt[,-c(127:129)], #Species untransformed abundances? vs transformed abundances?
                              expand = T) # The same weighted variance as the ordination scores.
spe.cmd.ne$species <- wascores(spe.cmd.ne$points, # Ordination scores for sites.
                              spe.dt[,-c(127:129)], #Species untransformed abundances? vs transformed abundances?
                              expand = T) # The same weighted variance as the ordination scores.
# Calculate species scores which are not calculated in cmdscale?
spe.h.cmd$species <- wascores(spe.h.cmd$points, # Ordination scores for sites.
                              spe.hel[,-c(127:129)], #Species untransformed abundances? vs transformed abundances?
                              expand = T) # The same weighted variance as the ordination scores.

par(mfrow = c(1,3))
ordiplot(spe.cmd, type = "text", main = c("Bray-Curtis Transformed Data"))
ordiplot(spe.h.cmd, type = "text", main = c("Hellinger Transformed Data"))
ordiplot(spe.cmd.ne, type = "text", main = c("Euclidean Data"))




#### Diversity ####
richness <- specnumber(x = spe.dt)
shan_div <- diversity(spe.dt, index = "shannon")
simp_div <- diversity(spe.dt, index = "simpson") #Calculates Gini-Simpson (1-D) 
invsimp_div <- diversity(spe.dt, index = "invsimpson") #
simp_evenness <- invsimp_div/richness #The effective number of species evenness


x <- as.data.frame(richness)
x$richness2 <- richness
par(mfrow = c(1,2))
plot(x$richness, invsimp_div)
abline(lm(x$richness ~ x$richness2), col = "red")
plot(richness, (invsimp_div/richness))


#Assign diversity measures 
spe.dt$richness <- richness; spe.dt$evenness <- simp_evenness
spe.dt$shan_div <- shan_div; spe.dt$simp_div <- simp_div
#Assign cmdscale outputs to species data set
spe.dt$cmd1 <- spe.cmd$points[,1]; spe.dt$cmd2 <- spe.cmd$points[,2]; spe.dt$cmd3 <- spe.cmd$points[,3]

par(mfrow = c(1,3))
eig_1 <- round(spe.cmd$eig*100/sum(spe.cmd$eig),2)[1]
eig_2 <- round(spe.cmd$eig*100/sum(spe.cmd$eig),2)[2]
eig_3 <- round(spe.cmd$eig*100/sum(spe.cmd$eig),2)[3]

plot(spe.cmd$points, type = "n", 
     main = "Species Beta Diversity - Community Composition",
     xlab = paste0("MDS-1 (", eig_1, "%)"),
     ylab = paste0("MDS-2 (", eig_2, "%)")
)
points(spe.cmd$points, cex = 2, pch = 21, col = "black", bg = 'darkred')
orditorp(spe.cmd, display = "species", 
         air = 1, cex = 1) #Surveyor labels are given space for readability 

plot(spe.cmd$points, type = "n", 
     main = "Species Beta Diversity - Community Composition",
     xlab = paste0("MDS-2 (", eig_2, "%)"),
     ylab = paste0("MDS-3 (", eig_3, "%)")
)
points(spe.cmd$points[,c(2,3)], cex = 2, pch = 21, col = "black", bg = 'darkred') 
orditorp(spe.cmd$species[,c(2,3)], display = "species", 
         air = 1, cex = 1) #Surveyor labels are given space for readability

plot(spe.cmd$points, type = "n", 
     main = "Species Beta Diversity - Community Composition",
     xlab = paste0("MDS-3 (", eig_3, "%)"),
     ylab = paste0("MDS-1 (", eig_1, "%)")
)
points(spe.cmd$points[,c(3,1)], cex = 2, pch = 21, col = "black", bg = 'darkred')
orditorp(spe.cmd$species[,c(3,1)], display = "species", 
         air = 1, cex = 1) #Surveyor labels are given space for readability

#Merge results with RRA_data
RRA_data <- merge(RRA_data, spe.dt, 
                  by.x = "RRA_id",
                  by.y = 0)

#### Landscape ####
#Match Landuse data set to new RRA_data
RRA_LU_filt %>% filter(RRA_id %in% 
         RRA_data$RRA_id) -> LU 

#Reclassify land uses to broad categories
LU %>% dplyr::select(Rural.Residential,
                              Estate.Residential, 
                              Medium.Density.Residential, 
                              High.Density.Residential) %>% #Maybe High.density residential should be urban...
  rowSums() -> RRA_data$residential
LU %>% dplyr::select(Industrial,
                              Commercial,
                              Institutional) %>% rowSums() -> RRA_data$urban
LU %>% dplyr::select(Golf.Course,
                              Recreational.Open.Space,
                              Agriculture, #Maybe Agriculture should be its own class? 
                              Vacant.Land) %>% rowSums() -> RRA_data$`green space`
LU %>% dplyr::select(Forest,
                              Successional.Forest,
                              Meadow,
                              Wetland,
                              Beach.Bluff) %>% 
  rowSums() -> RRA_data$natural
LU %>% dplyr::select(Roads,
                              Railway) %>% rowSums() -> RRA_data$transportation
LU %>% dplyr::select(Lacustrine,
                              Riverine) %>% rowSums() -> RRA_data$water

#Create Land Use PCA
#Will need to use a robust PCA? 

RRA_data %>% dplyr::select(urban, residential, transportation, 
                           `green space`, natural, water) -> LU_cat #Data for PCA
rownames(LU_cat) <- RRA_data$RRA_id #Assign rownames to unique id numbers of RRAs
RRA_LU_cor <- cor(LU_cat) #Determine correlation between land uses
RRA_LU_pca <- princomp(scale(LU_cat)) #Complete PCA with normalizing of data
#Extract pc1, pc2, pc3 for later analysis
RRA_data$pc1 <- RRA_LU_pca$scores[,1] 
RRA_data$pc2 <- RRA_LU_pca$scores[,2]
RRA_data$pc3 <- RRA_LU_pca$scores[,3]



#### Site Condition ####

## Habitat Types ##
#Convert project deliverable (desired ecotype)  to factor
RRA_data$Applicable.project.deliverable.categories. <- as.factor(RRA_data$Applicable.project.deliverable.categories.)

## Site level disturbances ##

#Invasive species
x = RRA_data$Invasive.species.observed.in.treated.location.
empty_list <- list()
for (i in 1:length(x)) {
  tmp <- stringr::str_split(x[i], pattern = ",")
  tmp <- sort(tmp[[1]])
  # print(tmp)
  empty_list[[i]] <- paste0(tmp[], collapse = " ")
}
RRA_data$invasive_species <- as.factor(unlist(empty_list))
#Invasive species richness
x = stringr::str_split(RRA_data$Invasive.species.observed.in.treated.location., ",")
y=lapply(x,function(l) {
  data.table(matrix(1,ncol=length(l),dimnames=list(1,l)))
})
y = rbindlist(y,fill=T)
setDF(y)
y[is.na(y)] = 0
y$`0`=NULL
rownames(y) = RRA_data$RRA_id
Invasive_species <- y
RRA_data$Obs_invSpecies <- rowSums(Invasive_species)

#Human interference factors
x <- RRA_data$Indicate.the.nature.of.human.interference.
empty_list <- list()
for (i in 1:length(x)) {
  tmp <- stringr::str_split(x[i], pattern = ",")
  tmp <- sort(tmp[[1]])
  empty_list[[i]] <- paste0(tmp[], collapse = " ")
}
RRA_data$human_interference <- as.factor(unlist(empty_list)) #Output is unique combinations of human disturbance factors
#Next we calculate the richness of human disturbances
x = stringr::str_split(RRA_data$Indicate.the.nature.of.human.interference., ",")
y=lapply(x,function(l) {
  data.table(matrix(1,ncol=length(l),dimnames=list(1,l)))
})
y = rbindlist(y,fill=T)
setDF(y)
y[is.na(y)] = 0
y$`0`=NULL
rownames(y) = RRA_data$RRA_id
human_interference <- y
RRA_data$NumHuman_Interf <- rowSums(human_interference) #Herbivore richness

#Herbivory factors
x <- RRA_data$Select.sources.of.browse.that.have.an.impact.on.project.success.
empty_list <- list()
for (i in 1:length(x)) {
  tmp <- stringr::str_split(x[i], pattern = ",")
  tmp <- sort(tmp[[1]])
  # print(tmp)
  empty_list[[i]] <- paste0(tmp[], collapse = " ")
}
RRA_data$browse_source <- as.factor(unlist(empty_list))
#Herbivory richness
x = stringr::str_split(RRA_data$Select.sources.of.browse.that.have.an.impact.on.project.success., ",")
y=lapply(x,function(l) {
  data.table(matrix(1,ncol=length(l),dimnames=list(1,l)))
})
y = rbindlist(y,fill=T)
setDF(y)
y[is.na(y)] = 0
y$`0`=NULL
rownames(y) = RRA_data$RRA_id
Browsing_source <- y
RRA_data$browsing_rich <- rowSums(Browsing_source) #Herbivore richness


#### Disturbance Management ####
### Need to figure out how to treat a time since treatment variable when there is a group with no treatment
RRA_data$Year.last.treated. <- RRA_data$Year.last.treated. #

#Human Exclusions
x <- RRA_data$Which.of.the.following.deterrants.are.present.
empty_list <- list()
for (i in 1:length(x)) {
  tmp <- stringr::str_split(x[i], pattern = " ")
  tmp <- sort(tmp[[1]])
  empty_list[[i]] <- paste0(tmp[], collapse = " ")
}
RRA_data$human_deterrants <- as.factor(unlist(empty_list)) 
# Human exclusions richness
x = stringr::str_split(RRA_data$Which.of.the.following.deterrants.are.present., ",")
y=lapply(x,function(l) {
  data.table(matrix(1,ncol=length(l),dimnames=list(1,l)))
})
y = rbindlist(y,fill=T)
setDF(y)
y[is.na(y)] = 0
y$`0`=NULL
rownames(y) = RRA_data$RRA_id
human_deterrents <- y
RRA_data$NumHuman_deterrents <- rowSums(human_deterrents)

#Herbivory deterrents
x <- fedmatch::clean_strings(RRA_data$Installations.present.to.deter.browse.)
x <- stringr::str_replace_all(x, c("n a" = "", "0" = "")); x <- trimws(x)
empty_list <- list()
for (i in 1:length(x)) {
  tmp <- stringr::str_split(x[i], pattern = " ")
  tmp <- sort(tmp[[1]])
  empty_list[[i]] <- paste0(tmp[], collapse = " ")
}
RRA_data$browse_deterrents <- as.factor(unlist(empty_list)) #creates vector of browsing deterrents
#Herbivory deterrents richness
x <- str_to_lower(RRA_data$Installations.present.to.deter.browse.)
x = stringr::str_split(x, ",")
y=lapply(x,function(x) {
  data.table(matrix(1,ncol=length(x),dimnames=list(1,x)))
})
y = rbindlist(y,fill=T)
setDF(y)
y[is.na(y)] = 0
#Remove unneeded columns
y$`0`=NULL; y$n_a=NULL
rownames(y) = RRA_data$RRA_id
#Set final dataset
Browsing_deterrents<- y
RRA_data$num_browsDeterrents <- rowSums(Browsing_deterrents)


#### Seasonality ####
RRA_data$planting_month <- month(RRA_data$planting_date_RRA) #What month the plantings were completed in?


#### Surveyor Bias ####
#Survey date (i.e. early vs late in the season and yearly differences?)
RRA_data$Survey_month <- ordered(month(RRA_data$Date.)) #What month the survey was completed?
RRA_data$Survey_year <- ordered(year(RRA_data$Date.)) #What year the survey was completed?


library(data.table)
x <- RRA_IBM$Survey_team #select survey_team information for analysis
x <- stringr::str_split(x, ",") #Separate by first single lower case with a space on either side and keep delimiter (?<=)
#Create data table of survyeors
y <- lapply(x, function(x) {
  data.table(matrix(1,ncol=length(x),dimnames=list(1,x)))
}
)
y = rbindlist(y,fill=T)
setDF(y)
y[is.na(y)] = 0
rownames(y) = RRA_IBM$rra_id
surv <- y #Save as survey dataset
library(vegan)
surv.dist <- vegdist(surv, method = "jaccard", binary = TRUE)

surv.cmd <- cmdscale(surv.dist, k = 2, 
                     eig = T, add = T)

surv.cmd$species <- wascores(surv.cmd$points,  # Ordination scores for sites.
                             surv, #Surveyors?
                             expand = T) # The same weighted variance as the ordination scores.
par(mfrow=c(1,1)); ordiplot(surv.cmd)

#Exploring NMDS as an option to capture surveyor bias
surv.nmds <- metaMDS(surv, 
                     distance = "jaccard", #Use Jaccard for Presence absence 
                     try = 20, trymax = 150, wascores = T)


surv.nmds
surv.nmds$stress

#Test of stress improvements
n = 5
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(surv, distance = "jaccard", k = i)$stress
}
names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)


par(mfrow=c(1,1))
plot(surv.nmds, type = 't', main =paste("NMDS/Bray - Stress =", round(surv.nmds$stress,3)),)

# devEMF::emf(file='NMDS_surveyors.emf',emfPlus = F)
par(mfrow=c(2,2))
#A - Stress plot across multiple dimensions
barplot(stress, ylab = "", 
        main = "Stress Test") 
abline(h = 0.01, col = "darkred")
#B #Shepard plot
spear <- round(cor(vegdist(surv, method = "jaccard"), 
                   dist(surv.nmds$points), 
                   method = "spearman"),3) #Spearman test
stressplot(surv.nmds, main="") 
text(x = 0.1, y =  6, #Hardcode the position of the spearman test text
     paste0("Spearman correlation = ", spear))
#C #paste("Goodness of Fit - NMDS/Jaccard - Stress =", round(surv.nmds$stress,3))
eig_1 <- round(surv.cmd$eig*100/sum(surv.cmd$eig),1)[1]
eig_2 <- round(surv.cmd$eig*100/sum(surv.cmd$eig),1)[2]
plot(surv.cmd$points, type = "n", 
     main = "",
     xlab = paste0("MDS-1 (", eig_1, "%)"),
     ylab = paste0("MDS-2 (", eig_2, "%)")
     )
points(surv.cmd$points, cex = 2, pch = 21, col = "black", bg = 'darkred') #Poorly fit sites have large bubbles
orditorp(surv.cmd, display = "species", 
         air = 1, cex = 1) #Surveyor labels are given space for readability 
#D
gof <- goodness(surv.nmds)
plot(surv.nmds$points, type = "n", 
     main = "", #paste("Goodness of Fit - NMDS/Jaccard - Stress =", round(surv.nmds$stress,3))
     xlab = paste0("NMDS-1"),
     ylab = paste0("NMDS-2")
     )
points(surv.nmds$points, pch = 21, col = "black", bg = 'steelblue',
         cex = gof*2/mean(gof) #Poorly fit sites have large bubbles
  )
orditorp(surv.nmds, display = "species", 
         air = 1, cex = 1 #Surveyor labels are given space for readability
         )
# dev.off()

# par(mfrow = c(1,1))
# ordiplot(surv.cmd, type = 'n') |>
#   points("sites", pch=21, col='red', bg='yellow') |>
#   text("species", col="blue", cex=0.9)
# 
# orditorp(surv.cmd, display = "sites", scaling = scl,
#          col = "blue", cex = 1, pch = 19)
# 
# df.surv <- as.data.frame(surv.cmd$points)
# 
# ggplot(df.surv, aes(x = V1, y = V2,
#        # color = "yellow",
#        # fill = "gray",
#        # palette = c('#1972A4', '#965F8A', '#FF7070', '#4AC6B7'),
#        add = "jitter"), ylab = "Dim 1", xlab = "Dim 2") +
#   geom_point(size=1) 
# 


################# Clustering by distance surveyors




################# End


#### PCoA For community data
#Extract eigenvectors from first axis for survyeor differences
# RRA_data$pcoa_surv1 <- surv.cmd$points[,1] 
# RRA_data$pcoa_surv2 <- surv.cmd$points[,2]
# RRA_data$pcoa_surv3 <- surv.cmd$points[,3]
RRA_data$surv_nmds1 <- surv.nmds$points[,1]
RRA_data$surv_nmds2 <- surv.nmds$points[,2]
#Add column for each surveyor for each abundant surveyor [1,0]
# RRA_data$Num_surveyors <- rowSums(y) #Create unique name
#Get the surveyors with the most differences in 
# surveyors <- colSums(y)
#Create common unique surveyors
RRA_data$surveyor1 <- y$`jennifer r`
RRA_data$surveyor3 <- y$`erica s`
RRA_data$surveyor5 <- y$`melissa s`
#Create surveyor_ids for most abundant survyeors and others
surveyorID <- RRA_data %>% dplyr::select(RRA_id, surveyor1, surveyor3, surveyor5) %>% 
  pivot_longer(cols = c(surveyor1:surveyor5), #Create new columns in long form for surveyor ids
               names_to = c("surveyor_id"),
               values_to = "surveyor") %>%
  filter(surveyor > 0) #remove rows without one of the 'big three' 
#Add new column back to original dataframe
RRA_data <- merge(RRA_data, surveyorID, by = "RRA_id", all.x = T)
RRA_data$surveyor_id[is.na(RRA_data$surveyor_id)] <- "surveyorO" #Reclassify NAs in surveyor_ids to surveyorO (others)




# Subset data -------------------------------------------------------------------
clmm_data <- RRA_data %>% dplyr::select(
  #Primary key
  RRA_id, 
  
  #Confounding Variables - Unmeasured Variability and groupings
  site_name, #Site history, unmeasured initial conditions.  
  Applicable.project.deliverable.categories., #Unaccounted for differences in project types
  # pcoa_surv1, pcoa_surv2, pcoa_surv3, #Team identity for restoration survey
  surv_nmds1, surv_nmds2,
  surveyor_id, #Unique surveyors
  Survey_team, #Unique teams
  
  #Response Variables
  ConifVegScore, DecidVegScore, ShrubVegScore, HerbVegScore, #Planted component scores
  Natural.cover.component.score., Natural.native.regeneration.score., #Site level metrics for NCC, NNR,   
  
  Survival.score., Survival.score..1, Survival.score..2, #Planted component survival scores
  
  Average.health.of.surviving.plants., Average.health.of.surviving.plants..1, #Planted component average health (leaf loss/browsing and discolouration)
  Average.health.of.surviving.plants..2, 
  
    
  ##Primary Predictors of Interest = Diversity metrics 
  richness, evenness, shan_div, simp_div,
  # pcoa1, pcoa2, pcoa3, 
  cmd1, cmd2, cmd3,
  
  #Predictors 0 - Temporal spatial
  Assessment.Year., #The number of years after planting (i.e., accumulation of mortality and stress).
  Survey_month, #Seasonality and accumulation of stress
  Survey_year, #Year of survey yearly climate and teams
  planting_month, #Seasonality of planting
  
  #Predictors 1 - Biotic disturbances 
  Browsing.of.planted.vegitation., Extent.of.browsing., Browse.Score., browse_source, browsing_rich, #Relative browsing Intensity metrics
  
  Invasive.species., Extent.of.invasive.species., Invasive.Score., Obs_invSpecies, #Relative invasive species abundances
  Human.interference.with.site., Extent.of.human.interference., Human.Score., human_interference, NumHuman_Interf,#Relative human distrubance 
  
  #Predictor 1b - Management
  Year.last.treated.,
  browse_deterrents, num_browsDeterrents,
  human_deterrants, NumHuman_deterrents,
  
  #Predictors 2 - Environmental/Landscape
  x, y, pc1, pc2, dist_TO, urban, residential, transportation, `green space`, natural, water,
  
  #Species abundances
  Abies.balsamea:Zizia.aurea
  
  
) 

#Scale numerics
clmm_data$dist_TO <- as.numeric(scale(clmm_data$dist_TO))
clmm_data$shan_div <- as.numeric(scale(clmm_data$shan_div))
clmm_data$simp_div <- as.numeric(scale(clmm_data$simp_div))
clmm_data$richness <- as.numeric(scale(clmm_data$richness))
clmm_data$evenness <- as.numeric(scale(clmm_data$evenness))
# clmm_data$pcoa1 <- as.numeric(scale(clmm_data$pcoa1))
# clmm_data$pcoa2<- as.numeric(scale(clmm_data$pcoa2))
clmm_data$Assessment.Year. <- as.numeric(as.character(clmm_data$Assessment.Year.))

clmm_data$urban <- as.numeric(scale(clmm_data$urban))
clmm_data$residential <- as.numeric(scale(clmm_data$residential))
clmm_data$transportation <- as.numeric(scale(clmm_data$transportation))
clmm_data$`green space` <- as.numeric(scale(clmm_data$`green space`))
clmm_data$natural <- as.numeric(scale(clmm_data$natural))
clmm_data$water <- as.numeric(scale(clmm_data$water))

#Disturbances
clmm_data$NumHuman_Interf <- ordered(clmm_data$NumHuman_Interf)
clmm_data$browsing_rich <- ordered(clmm_data$browsing_rich)
clmm_data$Obs_invSpecies <- ordered(clmm_data$Obs_invSpecies)
clmm_data$NumHuman_deterrents <- ordered(clmm_data$NumHuman_deterrents)
clmm_data$num_browsDeterrents <- ordered(clmm_data$num_browsDeterrents)
#Time since year treated is ordered with no treatment (0) worse then oldest treatment (2015),... and so on. 
clmm_data$Year.last.treated. <- ordered(as.numeric(clmm_data$Year.last.treated.))

#Create a summed disturbance measure 
#This may require new cut-points 
clmm_data$disturbance_intensity <- ordered(as.numeric(clmm_data$NumHuman_Interf) * 
                                             as.numeric(clmm_data$browsing_rich) * 
                                             as.numeric(clmm_data$Obs_invSpecies))
#Collapse Factors into new groups
clmm_data$disturbance_cat <- fct_collapse(clmm_data$disturbance_intensity, 
  disturbance_none = c("1"), #No disturbances observed at site.
  disturbance_v.low = c("2","3"), #One to two disturbances of a single type observed.
  disturbance_low = c("4"), #One to two disturbances of multiple types observed.
  disturbance_middle = c("5", "6", "8"), #Disturbances observed of multiple types
  disturbance_high = c("9", "10"), #Several disturbances observed of multiple types
  disturbance_v.high = c("12", "15", "16", "32") #Several or all disturbances observed of multiple types
)

# clmm_data$browse_deterrents
clmm_data$num_browsDeterrents
# clmm_data$human_deterrants 
clmm_data$NumHuman_deterrents
clmm_data$Year.last.treated.

#Identified sites that were treated more than two years prior to the survey.
#All treated sites were  within a five year window but about a third over two years old. 
inv_treat <- as.numeric(as.character(clmm_data$Survey_year)) - as.numeric(as.character(clmm_data$Year.last.treated.)) 
inv_treat 
#convert to binary treatment no treatment of invasive species in the last five years
clmm_data$inv_treatment <- fct_collapse(clmm_data$Year.last.treated.,
             untreated = c("0"),
             treated = c("2015", "2018", "2019", "2020", "2021", "2022")
)

clmm_data$deterent <- ordered(as.numeric(clmm_data$inv_treatment) * 
          as.numeric(clmm_data$num_browsDeterrents) * 
          as.numeric(clmm_data$NumHuman_deterrents))
clmm_data$deterent_cat <- fct_collapse(clmm_data$deterent,
                                       deterrents_none = c("1"), #No disturbance deterents
                                       deterrents_low = c("2", "3", "4"), #Some deterents for one - two disturbances
                                       deterrents_high = c("6", "8") #Multiple deterents for disturbances
                                       )
#Summary table of disturbances and deterrents
table(clmm_data$disturbance_cat, clmm_data$deterent_cat)

#Create working copy
clmm_data2 <- clmm_data
#Convert to characters to combine columns
clmm_data2$Average.health.of.surviving.plants. <- as.character(clmm_data2$Average.health.of.surviving.plants.)
clmm_data2$Average.health.of.surviving.plants..1 <- as.character(clmm_data2$Average.health.of.surviving.plants..1)
clmm_data2$Average.health.of.surviving.plants..2 <- as.character(clmm_data2$Average.health.of.surviving.plants..2)

## Survival ##
clmm_data2 %>% dplyr::rename(c("Coniferous" = "Survival.score.", 
                               "Deciduous" = "Survival.score..1",
                               "Shrub" = "Survival.score..2")) %>% 
  pivot_longer(cols = c(Coniferous:Shrub),
               names_to = c("planting_type"),
               values_to = "Survival_score") -> clmm_surv
clmm_surv$Survival_score <- as.ordered(clmm_surv$Survival_score) #Set to ordered factors
clmm_surv$planting_type <- as.factor(clmm_surv$planting_type)  #set new variable #planting type to factor

## Health ##
#Rename and combine columns into planting type and health score
clmm_data2 %>% dplyr::rename(c("Coniferous" = "Average.health.of.surviving.plants.", 
                               "Deciduous" = "Average.health.of.surviving.plants..1",
                               "Shrub" = "Average.health.of.surviving.plants..2")) %>% 
  pivot_longer(cols = c(Coniferous:Shrub),
               names_to = c("planting_type"),
               values_to = "Health_score") -> clmm_health
clmm_health$Health_score <- as.ordered(clmm_health$Health_score) #Convert to ordered factors
clmm_health$planting_type <- as.factor(clmm_health$planting_type) #set new variable planting type to factor


## Test of multi-collinarity ##
lm.survival<-lm(as.numeric(Survival_score) ~  hill_rich + hill_evenness + #Richness and evenness of plantings
                  cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                  pc1 + pc2 + dist_TO + #Landscape Factors
                  planting_month + #Dispersal/phenology 
                  # NumHuman_Interf + browsing_rich + Obs_invSpecies + #Disturbances
                  disturbance_cat + deterent_cat + 
                  Assessment.Year. + surv_nmds1 + surv_nmds2,
                data = clmm_surv[!is.na(clmm_surv$Survival_score),])

lm.health<-lm(as.numeric(Health_score) ~  hill_rich + hill_evenness + #Richness and evenness of plantings
                cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                pc1 + pc2 + dist_TO + #Landscape Factors
                planting_month + #Dispersal/phenology 
                # NumHuman_Interf + browsing_rich + Obs_invSpecies + #Disturbances
                disturbance_cat + deterent_cat +
                Assessment.Year. + surv_nmds1 + surv_nmds2,
              data = clmm_health[!is.na(clmm_health$Health_score),])

summary(lm.survival)
vif(lm.survival) 
summary(lm.health)
vif(lm.health)

# Set Labels for variables ----------------------------------------------------- 
library(sjlabelled)
## Survival Labels -------------------------------------------------------------
set_label(clmm_surv$planting_month) <- "Planting month"
set_label(clmm_surv$Assessment.Year.) <- "Assessment year"
set_label(clmm_surv$Obs_invSpecies) <-  "Invasive species richness"
set_label(clmm_surv$browsing_rich) <- "Herbivore richness"
set_label(clmm_surv$NumHuman_Interf) <- "Human disturbance"
set_label(clmm_surv$simp_div) <- "Simpsons diversity"
set_label(clmm_surv$shan_div) <- "Shannon diversity"
set_label(clmm_surv$richness) <- "Richness"
set_label(clmm_surv$evenness) <- "Evenness"
set_label(clmm_surv$cmd1) <- "Community composition PCo-1"
set_label(clmm_surv$cmd2) <- "Community composition PCo-2"
set_label(clmm_surv$cmd3) <- "Community composition PCo-3"
set_label(clmm_surv$pc1) <- "Urbanization PC-1"
set_label(clmm_surv$pc2) <- "Urbanization PC-2"
set_label(clmm_surv$dist_TO) <- "Distance to city centre"
set_label(clmm_surv$NumHuman_deterrents) <- "Human deterrents"
set_label(clmm_surv$num_browsDeterrents) <- "Herbivore deterrents"
set_label(clmm_surv$Year.last.treated.) <- "Invasive species treatments"
set_label(clmm_surv$surv_nmds1) <- "Surveyor bias 1"
set_label(clmm_surv$surv_nmds2) <- "Surveyor bias 2"
## Health Labels ---------------------------------------------------------------
set_label(clmm_health$planting_month) <- "Planting month"
set_label(clmm_health$Assessment.Year.) <- "Assessment year"
set_label(clmm_health$Obs_invSpecies) <-  "Invasive species richness"
set_label(clmm_health$browsing_rich) <- "Herbivore richness"
set_label(clmm_health$NumHuman_Interf) <- "Human disturbance"
set_label(clmm_health$simp_div) <- "Simpsons diversity"
set_label(clmm_health$shan_div) <- "Shannon diversity"
set_label(clmm_health$richness) <- "Richness"
set_label(clmm_health$evenness) <- "Evenness"
set_label(clmm_health$cmd1) <- "Community composition PCo-1"
set_label(clmm_health$cmd2) <- "Community composition PCo-2"
set_label(clmm_health$cmd3) <- "Community composition PCo-3"
set_label(clmm_health$pc1) <- "Urbanization PC-1"
set_label(clmm_health$pc2) <- "Urbanization PC-2"
set_label(clmm_health$dist_TO) <- "Distance to city centre"
set_label(clmm_health$NumHuman_deterrents) <- "Human deterrents"
set_label(clmm_health$num_browsDeterrents) <- "Herbivore deterrents"
set_label(clmm_health$Year.last.treated.) <- "Invasive species treatments"
set_label(clmm_health$surv_nmds1) <- "Surveyor bias 1"
set_label(clmm_health$surv_nmds2) <- "Surveyor bias 2"

## Data labels -----------------------------------------------------------------
set_label(clmm_data$planting_month) <- "Planting month"
set_label(clmm_data$Assessment.Year.) <- "Assessment year"
set_label(clmm_data$Obs_invSpecies) <-  "Invasive species richness"
set_label(clmm_data$browsing_rich) <- "Herbivore richness"
set_label(clmm_data$NumHuman_Interf) <- "Human disturbance"
set_label(clmm_data$simp_div) <- "Simpsons diversity"
set_label(clmm_data$shan_div) <- "Shannon diversity"
set_label(clmm_data$richness) <- "Richness"
set_label(clmm_data$evenness) <- "Evenness"
set_label(clmm_data$cmd1) <- "Community composition PCo-1"
set_label(clmm_data$cmd2) <- "Community composition PCo-2"
set_label(clmm_data$cmd3) <- "Community composition PCo-3"
set_label(clmm_data$pc1) <- "Urbanization PC-1"
set_label(clmm_data$pc2) <- "Urbanization PC-2"
set_label(clmm_data$dist_TO) <- "Distance to city centre"
set_label(clmm_data$NumHuman_deterrents) <- "Human deterrents"
set_label(clmm_data$num_browsDeterrents) <- "Number of Herbivore Deterrents"
set_label(clmm_data$browse_deterrents) <- "Herbivore Deterrents"
set_label(clmm_data$Year.last.treated.) <- "Invasive species treatments"
set_label(clmm_data$surv_nmds1) <- "Surveyor bias 1"
set_label(clmm_data$surv_nmds2) <- "Surveyor bias 2"

# Cumulative Link Mixed Models & Proportional Odds Regression Models -----------
library(ordinal)
library(emmeans)
## Natural Cover Component Analysis ---------------------------------------------- 
cm.NCC <- clmm(Natural.cover.component.score. ~  dist_TO + #distance to city center
      pc1 + pc2 + #Urbanization predict outcomes?
      disturbance_cat + #How do greater disturbances influence outcomes?
      richness + evenness + #Richness and evenness of plantings influence survival?
      cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival?
      planting_month + #Do the timing of plantings significantly influence restoration outcomes?
      deterent_cat + # and deterrents influence restoration outcomes?
      Assessment.Year. + #Do restorations improve in each year?
      surv_nmds1 + surv_nmds2 +
        (1|RRA_id),
    Hess = T,
    data = clmm_data); summary(cm.NCC)


## Survival Analysis -------------------------------------------------------------
cm.survival_data <- clmm_surv[!is.na(clmm_surv$Survival_score),]
cm.survival_data$RRA_id <- as.factor(cm.survival_data$RRA_id)
### Null model -------------------------------------------------------------------
cn.survival <- clmm(Survival_score ~
                      (1|Applicable.project.deliverable.categories.) + 
                      (1|Assessment.Year.), 
                    data = cm.survival_data); summary(cn.survival)


### Fixed effects model -----------------------------------------------------------
fix.survival <- clm(Survival_score ~  dist_TO + #distance to city center
                      pc1 + pc2 + #Urbanization predict outcomes?   
                      disturbance_cat + #How do greater disturbances influence outcomes?
                      richness + evenness + #Richness and evenness of plantings influence survival? 
                      cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                      planting_month + #Do the timing of plantings significantly influence restoration outcomes?  
                      deterent_cat + # and deterrents influence restoration outcomes? 
                      Assessment.Year. + #Do restorations improve in each year? 
                      surv_nmds1 + surv_nmds2,  #Survey
                    Hess = T,
                    data = cm.survival_data); summary(fix.survival)

### Mixed model ------------------------------------------------------------------
cm.survival <- clmm(Survival_score ~  pc1 + pc2 + #Urbanization predict outcomes?   
                      dist_TO + #distance to city center
                      disturbance_cat + #How do greater disturbances influence outcomes?
                      richness + evenness + #Richness and evenness of plantings influence survival? 
                      cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                      planting_month + #Do the timing of plantings significantly influence restoration outcomes?
                      deterent_cat + # and deterrents influence restoration outcomes? 
                      Assessment.Year. + #Do restorations improve in each year?
                      surv_nmds1 + surv_nmds2 + #Does differences in the survey team predict the reported restoration outcomes? 
                      (1|RRA_id) , #The survey id captures unmeasured elements of each survey such as 
                    Hess = T,
                    data = cm.survival_data); summary(cm.survival)


#### Test proportional odds assumption -----------------------------------------
nom.survival <- nominal_test(fix.survival) #Breach of proportional odds assumption.
write.csv(nom.survival, "data/1.raw/2.working/3.output/nominal_test_survival.csv")
# clmm2 which has approaches for specifying nominal effects where clmm does not
cm2.survival <- clmm2(Survival_score ~  pc1 + pc2 + #Urbanization predict outcomes?   
                        # dist_TO + #distance to city center
                        disturbance_cat + #How do greater disturbances influence outcomes?
                        richness + evenness + #Richness and evenness of plantings influence survival?
                        cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                        planting_month + #Do the timing of plantings significantly influence restoration outcomes?
                        deterent_cat + # and deterrents influence restoration outcomes? 
                        # Assessment.Year. + #Do restorations improve in each year?
                        surv_nmds1 + surv_nmds2, 
                      random = RRA_id, 
                      # nominal = ~ Assessment.Year.,
                      Hess = T,
                      data = cm.survival_data); summary(cm2.survival)
cm2n.survival <- clmm2(Survival_score ~  pc1 + pc2 + #Urbanization predict outcomes?   
                         # dist_TO + #distance to city center
                         disturbance_cat + #How do greater disturbances influence outcomes?
                         richness + evenness + #Richness and evenness of plantings influence survival?
                         cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                         planting_month + #Do the timing of plantings significantly influence restoration outcomes?
                         deterent_cat + # and deterrents influence restoration outcomes? 
                         Assessment.Year. + #Do restorations improve in each year?
                         surv_nmds1 + surv_nmds2, 
                       random = RRA_id, 
                       nominal = ~ Assessment.Year.,
                       Hess = T,
                       data = cm.survival_data); summary(cm2n.survival)
# Can test differences between models manually...


#### Reduced Model Selection ---------------------------------------------------
cm.survival.drop <- drop1(cm.survival)
cm.survival.drop[cm.survival.drop$AIC == min(cm.survival.drop$AIC), ]

f.survival <-  clm(Survival_score ~  pc2 + disturbance_cat +
                      evenness + cmd3 +
                      planting_month,
                   Hess = T,
                    data = cm.survival_data)
cmf.survival <- clmm(Survival_score ~  pc2 + disturbance_cat +
                       evenness + cmd3 +
                       planting_month +
                       (1|RRA_id), #Survey as random variable.
                     Hess = T,
                     data = cm.survival_data); summary(cmf.survival)

##### Test model differences ---------------------------------------------------
anova(fix.survival, cm.survival) #test mixed effect model differs (it does and is lower AIC) ***
anova(cmf.survival, cm.survival) # The reduced mixed effect model is not sig diff. but has lower AIC

# Health Analysis --------------------------------------------------------------
cm.health_data <- clmm_health[!is.na(clmm_health$Health_score),]
cm.health_data$RRA_id <- as.factor(cm.health_data$RRA_id)
## Null model ------------------------------------------------------------------
cn.health <- clmm(Health_score ~ (1|RRA_id),
                  data = cm.health_data)

## Fixed effects model ---------------------------------------------------------
fix.health <- clm(Health_score ~ dist_TO + #distance to city center
                    pc1 + pc2 + #Urbanization predict outcomes?   
                    disturbance_cat + #How do greater disturbances influence outcomes?
                    richness + evenness + #Richness and evenness of plantings influence survival? 
                    cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                    planting_month + #Do the timing of plantings significantly influence restoration outcomes?  
                    deterent_cat + # and deterrents influence restoration outcomes? 
                    Assessment.Year. + #Do restorations improve in each year? 
                    surv_nmds1 + surv_nmds2,  #Survey
                  Hess = T,
                  data = cm.health_data); summary(fix.health)

## Mixed effects model ---------------------------------------------------------
cm.health <- clmm(Health_score ~ pc1 + pc2 +  #Urbanization predict outcomes?   
                    dist_TO +
                    disturbance_cat + #How do greater disturbances influence outcomes?
                    richness + evenness + #Richness and evenness of plantings influence survival? 
                    cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                    planting_month + #Do the timing of plantings significantly influence restoration outcomes?  
                    deterent_cat + # and deterrents influence restoration outcomes? 
                    Assessment.Year. + #Do restorations improve in each year? 
                    surv_nmds1 + surv_nmds2 + #Does differences in the survey team predict the reported restoration outcomes? 
                    (1|RRA_id), #The survey id captures unmeasured elements of each survey such as
                  Hess = T,
                  data = cm.health_data); summary(cm.health)
#### Test proportional odds assumption -----------------------------------------
nom.health <- nominal_test(fix.health) #Breach of proportional odds assumption.
write.csv(nom.health, "data/1.raw/2.working/3.output/nominal_test_health.csv")

# clmm2 which has approaches for specifying nominal effects where clmm does not
cm2.health <- clmm2(Health_score ~  pc1 + pc2 + #Urbanization predict outcomes?   
                        dist_TO + #distance to city center
                        disturbance_cat + #How do greater disturbances influence outcomes?
                        richness + evenness + #Richness and evenness of plantings influence survival?
                        cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                        planting_month + #Do the timing of plantings significantly influence restoration outcomes?
                        deterent_cat + # and deterrents influence restoration outcomes? 
                        # Assessment.Year. + #Do restorations improve in each year?
                        surv_nmds1 + surv_nmds2, 
                      random = RRA_id, 
                      # nominal = ~ Assessment.Year.,
                      Hess = T,
                      data = cm.health_data); summary(cm2.health)
cm2n.health <- clmm2(Health_score ~  pc1 + pc2 + #Urbanization predict outcomes?   
                         dist_TO + #distance to city center
                         disturbance_cat + #How do greater disturbances influence outcomes?
                         richness + evenness + #Richness and evenness of plantings influence survival?
                         cmd1 + cmd2 + cmd3 + #Does what is planted influence the outcome of the survival? 
                         planting_month + #Do the timing of plantings significantly influence restoration outcomes?
                         deterent_cat + # and deterrents influence restoration outcomes? 
                         Assessment.Year. + #Do restorations improve in each year?
                         surv_nmds1 + surv_nmds2, 
                       random = RRA_id, 
                       nominal = ~ Assessment.Year.,
                       Hess = T,
                       data = cm.health_data); summary(cm2n.health)
# Can test differences between models manually...

#### Reduced Model Selection ---------------------------------------------------
cm.survival.drop <- drop1(cm.survival)
cm.survival.drop[cm.survival.drop$AIC == min(cm.survival.drop$AIC), ]

f.survival <-  clm(Survival_score ~  pc2 + disturbance_cat +
                     evenness + cmd3 +
                     planting_month,
                   Hess = T,
                   data = cm.survival_data)
cmf.survival <- clmm(Survival_score ~  pc2 + disturbance_cat +
                       evenness + cmd3 +
                       planting_month +
                       (1|RRA_id), #Survey as random variable.
                     Hess = T,
                     data = cm.survival_data); summary(cmf.survival)

##### Test model differences ---------------------------------------------------
anova(fix.survival, cm.survival) #test mixed effect model differs (it does and is lower AIC) ***
anova(cmf.survival, cm.survival) # The reduced mixed effect model is not sig diff. but has lower AIC


#### Reduced Model Selection ---------------------------------------------------
cm.health.drop <- drop1(cm.health, trace = T)
cm.health.drop[cm.health.drop$AIC == min(cm.health.drop$AIC), ]

f.survival <-  clm(Survival_score ~  pc2 + disturbance_cat +
                     evenness + cmd3 +
                     planting_month,
                   Hess = T,
                   data = cm.survival_data)
cmf.health <- clmm(Health_score ~ cmd3 + surv_nmds2 +
                        (1|RRA_id), #The survey id captures unmeasured elements of each survey such as 
                      Hess = T,
                      data = clmm_health[!is.na(clmm_health$Health_score),]); summary(cmf.health)

anova(cm.health, cn.health) #Fixed effects **
anova(cm.health, fix.health) #Random effects ***
anova(cm.health.min, cn.health)
stepf_health <- step(fix.health, direction = "both")

# Minimal model
min.health <- clmm(Health_score ~ hill_shan + cmd3 + pc1 + pc2 + planting_month + 
                     NumHuman_Interf + browsing_rich + Obs_invSpecies + planting_type + 
                     pcoa_surv1 +
                       (1|RRA_id), #Site as random variable.
                     Hess = T,
                     data = clmm_health[!is.na(clmm_health$Health_score),])
summary(min.health)



#### Drop1 - Likelihood Ratio Tests #### 




fix.survival.step <- step(fix.survival)
fix.health.step <- step(fix.health)
cm.survival.step <- clmm(Survival_score ~ pc1 + pc2 + disturbance_cat + 
                           hill_evenness + cmd3 + 
                           planting_month + 
                           (1|RRA_id), 
                         Hess = T,
                         data = cm.survival_data); summary(cm.survival.step)
cm.health.step <- clmm(Health_score ~ pc1 + pc2 + disturbance_cat +
                         cmd3 + planting_month + 
                         deterent_cat + Assessment.Year. + 
                         surv_nmds1 + 
                           (1|RRA_id), 
                         Hess = T,
                         data = cm.health_data); summary(cm.health.step)
###
sink("drop1_cmsurvival.txt") #Export the r output to text file
drop1(cm.survival, test = "Chi")
sink()
cm.survival.drop <- clmm(Survival_score ~ pc2 +
                           hill_evenness + cmd3 + 
                           planting_month + 
                           (1|RRA_id), 
                         Hess = T,
                         data = cm.survival_data)
anova(cm.survival, cm.survival.drop)

drop1(cm.health, test = "Chi")
cm.health.drop <- clmm(Health_score ~ pc1 + pc2 + disturbance_cat + 
                         planting_month + 
                         deterent_cat + 
                         Assessment.Year. + 
                         surv_nmds1 + 
                        (1|RRA_id), 
                      Hess = T,
                      data = cm.health_data); summary(cm.health.drop)
anova(cm.health, cm.health.drop)


anova(cm.survival.step, cm.survival)
anova(cm.health.step, cm.health)
#### Figures ####


### Main Text ###
library(sjPlot)

## Results ##
par(mfrow = c(1,1))
## Land use variation ##
biplot(RRA_LU_pca, 
       data = LU_cat, 
       colour = 'RRA_id',
       col=c('lightblue', 'darkred'),
       cex=c(1, 1.3),
       xlim=c(-.15, .15),
       main='',
       xlab='First Component (27.7%)',
       ylab='Second Component (25.7)',
       expand=1.2)
# devEMF::emf(file='Figure_Landuse_PCA.emf',emfPlus = F)
# LU_biplot
# dev.off()

library(vegan3d)
ordirgl(spe.cmd$species[,1:3], type = 't')

##### Come back here to fix plots #####
# Fitness #
#Survival
par(mfrow = c(1,2))
# devEMF::emf(file='figure/Figure_survival.emf',emfPlus = F)
plot_model(cm.survival, type = "est", 
           rm.terms = c("0|1", "1|2", "2|3", "3|4"), #Remove intercepts
           vline.color = "blue",
           show.values=TRUE,
           value.offset = .3,
           show.p=TRUE,
           line.size = 1,
           dot.size = 3,
           title="",
           show.data = T)
# dev.off()
#
devEMF::emf(file='figure/Figure_health.emf',emfPlus = F)
plot_model(cm.health, type = "est", 
           rm.terms = c("0|1", "1|2", "2|3", "3|4"), #Remove intercepts
           vline.color = "blue",
           show.values=TRUE,
           value.offset = .3,
           show.p=TRUE,
           line.size = 1,
           dot.size = 3,
           title="",
           show.data = T) 
dev.off()

tab_model(cm.survival, cm.health,
          digits = 6,
          dv.labels = c("Survival","Health"),
          show.intercept = F,
          show.est = T
          # CSS = list(
          #   css.depvarhead = 'color: purple;',
          #   css.centeralign = 'text-align: left;', 
          #   css.firsttablecol = 'font-weight: bold;', 
          #   css.summary = 'color: blue;'),
          # file = "Survival-Health_Full.html"
          )
library(webshot)
webshot("Survival-Health_Full.html", "tab_model.png")

survival_plot <- plot_model(cm.survival,
                            type = "est",
                            terms = c("hill_rich","hill_evenness","cmd1","cmd2",           
                                      "cmd3", "pc1","pc2","dist_TO",
                                      "planting_month","disturbance_cat","deterent_cat","surv_nmds1",     
                                      "surv_nmds2"),
                            # axis.labels = c("Species Richness",
                            #                 "Species Evenness",
                            #                 "Community composition PCo-1","Community composition PCo-2","Community composition PCo-3",
                            #                 "Urbanization PC-1","Urbanization PC-2",
                            #                 "Distance to city-center",
                            #                 "Planting month",
                            #                 "Undisturbed","Little disturbance","Somewhat disturbed","Moderately disturbed","Disturbed","Heavily disturbed", 
                            #                 "No deterrents","Some deterrents","Deterrents",
                            #                 "Surveyor Bias NMDS-1","Surveyor Bias NMDS-2"),
                            vline.color = "blue",
                            show.intercept = F,
                            show.values=TRUE,
                            value.offset = .3,
                            show.p=TRUE,
                            line.size = 1,
                            dot.size = 3,
                            title="",
                            show.data = T)
devEMF::emf(file='Figure_survival.emf',emfPlus = F)
survival_plot
dev.off()
#Health
health_plot <- plot_model(min.health,
                          # group.terms = c(1,1,1,1,2,2,2,3,4,4,4,3,3),
                          terms = c("hill_shan","cmd3",
                                    "pc1","pc2",
                                    "planting_month", 
                                    "NumHuman_Interf","browsing_rich","Obs_invSpecies", 
                                    "Assessment.Year.", "planting_typeShrub", "planting_typeDeciduous", "pcoa_surv1"),
                          auto.label = T,
                          vline.color = "blue",
                          show.intercept = F,
                          show.values=TRUE,
                          value.offset = .3,
                          show.p=TRUE,
                          line.size = 1,
                          dot.size = 3,
                          title="",
                          show.data = T)
devEMF::emf(file='Figure_health.emf',emfPlus = F)
health_plot
dev.off()



c.survival <- clm(Survival_score ~ 1, 
                  data = cm.survival_data); summary(c.survival)

anova(cm.survival, cn.survival, type = "II")
anova(cm.health, cn.health, type = "II")

anova(cm.health.min, cn.health, type = "II")

sjPlot::tab_model(cm.survival, cm.health,
                  digits = 6,
                  #    pred.labels =c("(Intercept) - 10% Survival", "(Intercept) - 25% Survival", 
                  # "(Intercept) - 50% Survival", "(Intercept) - 75% Survival",
                  # "Simpson diversity", "Community composition Axis #1", 
                  # "Community composition Axis #2", "Community composition Axis #3",
                  # "Urbanization (PC1)", "Impervious land uses (PC2)", "Distance to city centre",
                  # "Planting month",
                  # "Human interferrence richness", "Herbivore richness", "Invasive species richness", 
                  # "Assessment year", "Surveyor bias"),
                  dv.labels= c("Survival", "Health"), 
                  show.intercept = F,
                  show.df = F, #Inf? why
                  
)




## Surveyor composition Figures ##

par(mfrow = c(1,1))
biplot.pcoa(surv.cmd, surv, 
            plot.axes = c(1,2),
            main = "", 
            col=c('lightblue', 'darkred'),
            cex=c(1, 1.3),
            xlim=c(-1, 1),
            # xlab='First Component (27.1%)',
            # ylab='Second Component (26.3)',
            expand=1.2) #This plot requires improvement





######## Supplementary #################################################
library(glue)
#S1
# #### For Visualizing ordination ##
# positions <- spe.h.cmd$points
# colnames(positions) <- c("PCo1", "PCo2", "PCo3")
# perc_expl <- 100 *spe.h.cmd$eig/sum(spe.h.cmd$eig)
# pretty_pe <- format(round(perc_expl <- 100 *spe.h.cmd$eig/sum(spe.h.cmd$eig), digits = 2), 
#                     nsmall=2,
#                     trim = T)
#Differences between hellinger transformed species composition ordination and untransformed
par(mfrow = c(1,3))
ordiplot(spe.h.cmd, type = "text", main = c("Hellinger Transformed Data"))
ordiplot(spe.cmd, type = "text", main = c("Bray-Curtis Transformed Data"))
ordiplot(spe.h.cmd2, type = "text", main = c("Hellinger Bray-Curtis Data"))

#Data set for clusters and functional types
spe.func <- read.csv("Data/species_functypes.csv")
#rename column to match for following merge
spe.func <- spe.func %>% rename(species_code = Species.Code)
########Grouping##############
spe_grouping <- merge(spe.func, spe.cmd$species, 
                      by.x = "species_code", by.y = 0, #The number zero specifies rownames?
                      all.y = T) #Keep all rows including none matches
library(naniar) #For NA replacement
spe_grouping <- spe_grouping %>%
  replace_with_na_all(condition = ~.x %in% common_na_strings)
spe_grouping <- spe_grouping %>% mutate(Sun.Exposure = tolower(Sun.Exposure)) %>% #Lowercase all entries
  mutate(Sun.Exposure = str_squish(Sun.Exposure)) %>% #Remove extra white spaces around strings and double spaces between characters
  mutate(Sun.Exposure = str_replace_all(Sun.Exposure," - ", "-")) #Remove depulicates due to space around " - "
spe_grouping[,4:9] <- lapply(spe_grouping[,4:9], function(x) {factor(x) #Set to factors 
}
)

#calculate the distances between species in the ordination
mds_spe.dist <- dist(na.omit(spe.cmd$species))
mds_spe.clust <- hclust(mds_spe.dist, method = "ward.D2")


#Create functional groups
# spe.func <- data.frame(colnames(spe.dt[,-c(1,257:262)]))
# spe.func <- rename(spe.func, "species_short" = "colnames.spe.dt....c.1..257.262...") 
# write.csv(spe.func, file = "species_functionalTypes.csv")

# biplot.pcoa(spe.h.pcoa, spe.dt[,-c(1, 257:259)],
#             xlim = c(-1,1),
#             plot.axes = c(1,3))

ordiplot(spe.h.cmd)


positions <- spe.h.cmd$points
colnames(positions) <- c("PCo1", "PCo2", "PCo3")
perc_expl <- 100 *spe.cmd$eig/sum(spe.cmd$eig)
pretty_pe <- format(round(perc_expl <- 100 *spe.h.cmd$eig/sum(spe.h.cmd$eig), digits = 2), 
                    nsmall=2,
                    trim = T)
cmd1_x <- spe.cmd$points[,1]
cmd2_y <- spe.cmd$points[,2]
cmd3_z <- spe.cmd$points[,3]

plot(x = cmd1_x, y = cmd2_y)

  geom_point() +
  theme_bw()

labs <- c(glue("PCo-1 ({pretty_pe[1]}%)"),
          glue("PCo-2 ({pretty_pe[2]}%)"))
positions <- spe.h.pcoa$vectors
positions %>% 
  as_tibble() %>% 
  ggplot(aes(x = Axis.1, y = Axis.2)) + 
  geom_point() +
  theme_bw() +
  labs(x = labs[1], y = labs[2])

tibble(pe = perc_expl, 
       axis = 1:length(perc_expl)) %>% 
  ggplot(aes(x = axis, y = pe)) + 
  geom_line() + 
  coord_cartesian(xlim = c(1,10))
fit.spe <- envfit(ord=spe.h.cmd, env=spe.dt[,-c(1, 257:259)])
plot(fit.spe, p.max = 0.05)
                           


biplot(spe.h.pcoa, spe.dt[,-c(1, 257:259)],
       plot.axes = c(1,2),
       
       main = "", 
       cex=c(5, 5.3),
       
       # xlim=c(-1, 1),
       expand= 2) #This plot requires improvement

biplot(spe.h.cmd, spe.dt[,-c(1, 257:259)],
       plot.axes = c(1,3),
       
       main = "", 
       cex=c(5, 5.3),
       
       # xlim=c(-1, 1),
       expand= 2) #This plot requires improvement


#S



library(ggfortify)

positions <- surv.cmd$points
colnames(positions) <- c("PCo1", "PCo2", "PCo3")
perc_expl <- 100 *surv.cmd$eig/sum(surv.cmd$eig)
pretty_pe <- format(round(perc_expl <- 100 *surv.cmd$eig/sum(surv.cmd$eig), digits = 2), 
                    nsmall=2,
                    trim = T)




ggplot(data = data.frame(positions),
       aes(x = PCo1, y = PCo2)) +
  geom_point() +
  theme_bw()
labs <- c(glue("PCo-1 ({pretty_pe[1]}%)"),
          glue("PCo-2 ({pretty_pe[2]}%)"))
positions %>% 
  as_tibble() %>% 
  ggplot(aes(x = PCo1, y = PCo2)) + 
  geom_point() +
  theme_bw() +
  labs(x = labs[1], y = labs[2])



labs <- c(glue("PCo-1 ({pretty_pe[1]}%)"),
          glue("PCo-2 ({pretty_pe[2]}%)"))
positions <- spe.h.pcoa$vectors
positions %>% 
  as_tibble() %>% 
  ggplot(aes(x = Axis.1, y = Axis.2)) + 
  geom_point() +
  theme_bw() +
  labs(x = labs[1], y = labs[2])


plot(spe.dt$cmd1, spe.dt$cmd3, type = "n", 
     axes = TRUE, xlab = "PCoA #1", ylab = "PCoA #3", 
     main = "")
text(spe.dt$cmd1, spe.dt$cmd3, rownames(spe.dt), cex = 0.6)

pl <- ordiplot(spe.h.cmd, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")

text(pl, "species", arrows = T, col="blue", cex=0.9)

p2 <- orditorp(spe.h.cmd, display = "species", scaling = scl,
               col = "forestgreen", pch = 2, cex = 1)

pl <- ordiplot(spe.h.cmd, type = "none")

ordilabel(spe.h.cmd, display = "species", groups = spe_grouping$Functional.Group)

points(pl, "species", pch=19, col=spe_grouping$Functional.Group)

ordipointlabel(spe.h.cmd, scaling = 3)
ordipointlabel(spe.h.cmd, display = "species", scaling = "symm", 
               add = TRUE)

######################



spe_grouping <- merge(spe.func, spe.h.cmd$species, 
                      by.x = "species_code", by.y = 0, #The number zero specifies rownames?
                      all.y = T) #Keep all rows including none matches
library(naniar) #For NA replacement
spe_grouping <- spe_grouping %>%
  replace_with_na_all(condition = ~.x %in% common_na_strings)
spe_grouping <- spe_grouping %>% mutate(Sun.Exposure = tolower(Sun.Exposure)) %>% #Lowercase all entries
  mutate(Sun.Exposure = str_squish(Sun.Exposure)) %>% #Remove extra white spaces around strings and double spaces between characters
  mutate(Sun.Exposure = str_replace_all(Sun.Exposure," - ", "-")) #Remove depulicates due to space around " - "
spe_grouping[,4:9] <- lapply(spe_grouping[,4:9], function(x) {factor(x) #Set to factors 
}
)



# spe.h.cmd$species <- merge(spe.h.cmd$species, spe_grouping[,1:9], 
#                            by.x = 0, by.y = "species_code", #The number zero specifies rownames?
#                            all = F) #Keep all rows including none matches


# #Filling missing data on species planted regarding functional groups, genus, etc.
# Genus <- spe.dt
# tmp_matching <- spe_grouping %>% select(species_code, Functional.Group)

########################33
# #Rename colnames to reduced detail functional groups
# names(Genus) <- tmp_matching$Functional.Group[match(names(Genus),tmp_matching$species_code)]
# ##combine columns by summing across 
# #*****Find code from cleaning planting lists for summing across multiple matching columns.
# # sapply(unique(colnames(Genus)), function(x) rowSums(Genus[,grepl(, colnames(Genus))]))
# Genus %>% group_by(colnames())

#Calculate clusters of species plantings
mds_spe.dist <- dist(na.omit(spe.cmd$species),
                     method = "euclidean") #Default euclidean method selected but need to determine its suitability
mds_spe.clust <- hclust(mds_spe.dist, # Distance matrix of WA scores
                        method = "ward.D2") # Need to determine correctness/suitability of this method vs alternates
# Determine cutpoints for clusters
library(cluster)
library(factoextra)
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(spe.cmd$species, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

# devEMF::emf(file='kmeans_communitycomposition.emf',emfPlus = F)
plot(k.values, wss_values,
     type="b", 
     pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
)
abline(v=6, col="blue", lty=c(2))
abline(v=8, col="darkred", lty=c(3))

# dev.off()
gap_kmeans <- clusGap(na.omit(spe.cmd$species), #Calculates a gap statistic to determine the optimal number of clusters 
                      FUN = kmeans, #Function to use is kmeans. What other options? kmeans/pam/
                      nstart = 30, #How many random sets should be chosen? What does this mean? 
                      # iter.max = 10, #Max number of iterations. What does this do?
                      K.max = 10, #Number of clusters to compute .. with greater than 10 we see a new local max calculated.
                      B = 100) #Number of Monte Carlo bootstrap sample (higher better but longer time?)
plot(gap_kmeans, main = "clusGap(., FUN = kmeans, n.start= 30, B= 100)")

# library(fpc)
# spe_dbscan <- dbscan(data = mds_spe.dist,
#                      eps = 0.2)
# spe_dbscan

#More robust version of k-means but much slower
final <- pam(x = na.omit(spe.cmd$species), k = 6)
gap_pam <- clusGap(na.omit(spe.cmd$species), #Calculates a gap statistic to determine the optimal number of clusters 
                   FUN = pam, #Function to use is kmeans. What other options? kmeans/pam/
                   nstart = 30, #How many random sets should be chosen? What does this mean? 
                   K.max = 10, #Number of clusters to compute .. with greater than 10 we see a new local max calculated.
                   B = 100) #Number of Monte Carlo bootstrap sample (higher better but longer time?)
plot(gap_pam, main = "clusGap(., FUN = kmeans, n.start= 30, B= 100)")



# devEMF::emf(file='gapstat_communitycomposition.emf',emfPlus = F)
fviz_nbclust(x = spe.cmd$species,
             FUNcluster = cluster::pam,
             method = "wss",
             k.max = 10, 
             nboot = 100 #Default value
             )
fviz_gap_stat(gap_pam, 
              maxSE = list(method = "globalmax") # Justification globalmax? vs some SEmax
) 
# dev.off()

# final <- kmeans(mds_spe.dist, #Numeric dissimilarity matrix
#                 centers = 8, # The number of clusters to be computed
#                 nstart = 30) #How many random sets should be chosen? What does this mean? 
#                 
library(viridis)

# devEMF::emf(file = "clusterplot_communitycomposition_MDS1-2.emf", emfPlus = F)
fc_12 <- fviz_cluster(final, data = mds_spe.dist,
                      # palette = c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951', 'black'),
                      axes = c(1,2),
                      repel =TRUE,
                      ellipse.type = "confidence",
                      main = "",
                      ggtheme =theme_minimal()
                      )  + scale_shape_manual('Cluster', values=c(20,21,22,23,24,25)
                      ) + scale_color_manual(values = palette.colors(n = 6, palette = "Dark 2")) 
# dev.off()
# devEMF::emf(file = "clusterplot_communitycomposition_MDS1-3.emf", emfPlus = F)
fc_13 <- fviz_cluster(final, data = mds_spe.dist,
                      axes = c(1,3),
                      repel =TRUE,
                      ellipse.type = "confidence",
                      main = "",
                      ggtheme =theme_minimal()
                      ) + scale_shape_manual('Cluster', values=c(20,21,22,23,24,25)
                      ) + scale_color_manual(values = palette.colors(n = 6, palette = "Dark 2"))  
# dev.off()
# devEMF::emf(file = "clusterplot_communitycomposition_MDS2-3.emf", emfPlus = F)
fc_23 <- fviz_cluster(final, data = mds_spe.dist,
                      axes = c(2,3),
                      repel =TRUE,
                      ellipse.type = "confidence",
                      main = "",
                      ggtheme =theme_minimal()
                      ) + scale_shape_manual('Cluster', values=c(20,21,22,23,24,25)
                      ) + scale_color_manual(values = palette.colors(n = 6, palette = "Dark 2"))  
# dev.off()
fc_12
fc_13
fc_23



spe.dt2 <- spe.dt %>% rownames_to_column('RRA_id')
RRA_data <- merge(RRA_data, spe.dt2[c(1,126:130)], by = "RRA_id")

######################
spe.func$species_code <- tolower(spe.func$species_code)
rownames(spe.cmd$species) <- tolower(rownames(spe.cmd$species))

spe_grouping <- merge(spe.func, spe.cmd$species, 
                      by.x = "species_code", by.y = 0, #The number zero specifies rownames?
                      all.y = T) #Keep all rows including none matches
library(naniar) #For NA replacement
spe_grouping <- spe_grouping %>%
  replace_with_na_all(condition = ~.x %in% common_na_strings)
spe_grouping <- spe_grouping %>% mutate(Sun.Exposure = tolower(Sun.Exposure)) %>% #Lowercase all entries
  mutate(Sun.Exposure = str_squish(Sun.Exposure)) %>% #Remove extra white spaces around strings and double spaces between characters
  mutate(Sun.Exposure = str_replace_all(Sun.Exposure," - ", "-")) #Remove depulicates due to space around " - "

library(plyr)
spe_grouping$Functional.Group <- plyr::revalue(spe_grouping$Functional.Group, c("Nt Tree" = "deciduous tree", 
                                                                                "Nt P-Grass" = "herbaceous plant",
                                                                                "Nt P-Forb" = "herbaceous plant",
                                                                                "Nt P-Sedge" = "herbaceous plant",
                                                                                "Nt Shrub" = "deciduous shrub",
                                                                                "Nt W-Vine" = "deciduous shrub") #Riverbank Grape 
)
spe_grouping$Coefficient.of.Wetness <- as.numeric(spe_grouping$Coefficient.of.Wetness)
spe_grouping[,c(3:6,8,9)] <- lapply(spe_grouping[,c(3:6,8,9)], function(x) {factor(x) #Set to factors 
}
)


# spe_grouping <- spe_grouping[-duplicated(spe_grouping$species_code),] #Remove duplicated row
#add clusters to species functional groups table
spe_grouping$clusters <- final$clustering

t1<- table(spe_grouping$Coefficient.of.Wetness, spe_grouping$clusters)
t1

spe_grouping$CoW <- as.numeric(scale(as.numeric(as.character(spe_grouping$Coefficient.of.Wetness))))
# aov.clus <- aov(as.factor(clusters) ~ Functional.Group, data = spe_grouping)
# summary(lm.clus) #

#cmd scores to CoW
V1.lm <- lm(V1 ~ CoW + Functional.Group, data = spe_grouping); summary(V1.lm)
V2.lm <- lm(V2 ~ CoW + Functional.Group, data = spe_grouping); summary(V2.lm)
V3.lm <- lm(V3 ~ CoW + Functional.Group, data = spe_grouping); summary(V3.lm)

library(car)
V1.anova <- Anova(V1.lm); V1.anova
V2.anova <- Anova(V2.lm); V2.anova
V3.anova <- Anova(V3.lm); V3.anova

# TukeyHSD(clus.aov, conf.level = .95)

library(emmeans)
#pairwise comparison 
V1.emm <- emmeans(V1.lm, ~Functional.Group)
V1.emm.df <- as.data.frame(V1.emm)
V1.emm.p <- pairs(emmeans(V1.lm, ~Functional.Group)) ### Use this for plotting differences between groups
v1.emm.p.df <- as.data.frame(V1.emm.p)

V2.emm <- emmeans(V2.lm, ~Functional.Group)
V2.emm.df <- as.data.frame(V2.emm)
V2.emm.p <- pairs(emmeans(V2.lm, ~Functional.Group)) ### Use this for plotting differences between groups
v2.emm.p.df <- as.data.frame(V2.emm.p)

V3.emm <- emmeans(V3.lm, ~Functional.Group)
V3.emm.df <- as.data.frame(V3.emm)
V3.emm.p <- pairs(emmeans(V3.lm, ~Functional.Group)) ### Use this for plotting differences between groups
v3.emm.p.df <- as.data.frame(V3.emm.p)


library(tidyr)
library(rstatix)
#V1
bxp_data <- separate_wider_delim(v1.emm.p.df, 
                                 cols = contrast, delim = " - ", 
                                 names = c("group1", "group2"))
bxp_data[bxp_data$p.value<=0.001,"p.adj.signif"]<-"***"
bxp_data[bxp_data$p.value>0.001 & bxp_data$p.value<=0.01,"p.adj.signif"]<-"**"
bxp_data[bxp_data$p.value>0.01 & bxp_data$p.value<=0.05,"p.adj.signif"]<-"*"
bxp_data[bxp_data$p.value>0.05 & bxp_data$p.value<=0.1, "p.adj.signif"]<-"."
bxp_data[bxp_data$p.value>0.1,"p.adj.signif"]<-"ns"
bxp_data.V1 <- bxp_data
bxp_data.V1$dimension <- "V1"
#V2
bxp_data <- separate_wider_delim(v2.emm.p.df, 
                                 cols = contrast, delim = " - ", 
                                 names = c("group1", "group2"))
bxp_data[bxp_data$p.value<=0.001,"p.adj.signif"]<-"***"
bxp_data[bxp_data$p.value>0.001 & bxp_data$p.value<=0.01,"p.adj.signif"]<-"**"
bxp_data[bxp_data$p.value>0.01 & bxp_data$p.value<=0.05,"p.adj.signif"]<-"*"
bxp_data[bxp_data$p.value>0.05 & bxp_data$p.value<=0.1, "p.adj.signif"]<-"."
bxp_data[bxp_data$p.value>0.1,"p.adj.signif"]<-"ns"
bxp_data.V2 <- bxp_data
bxp_data.V2$dimension <- "V2"
#V3
bxp_data <- separate_wider_delim(v3.emm.p.df, 
                                 cols = contrast, delim = " - ", 
                                 names = c("group1", "group2"))
bxp_data[bxp_data$p.value<=0.001,"p.adj.signif"]<-"***"
bxp_data[bxp_data$p.value>0.001 & bxp_data$p.value<=0.01,"p.adj.signif"]<-"**"
bxp_data[bxp_data$p.value>0.01 & bxp_data$p.value<=0.05,"p.adj.signif"]<-"*"
bxp_data[bxp_data$p.value>0.05 & bxp_data$p.value<=0.1, "p.adj.signif"]<-"."
bxp_data[bxp_data$p.value>0.1,"p.adj.signif"]<-"ns"
bxp_data.V3 <- bxp_data
bxp_data.V3$dimension <- "V3"

bxp_data <- rbind(bxp_data.V1, bxp_data.V2, bxp_data.V3)

bxp_spe_grouping <- spe_grouping %>% pivot_longer(cols = c(V1:V3), 
                                                  names_to = "dimension",
                                                  values_to = "eigenvalue")

bxp_df <- merge(bxp_data, bxp_spe_grouping[,c(1,6,11,12)], 
                by = "dimension") 

#Emmeans post-hoc test for functional groups
# ggboxplot(bxp_spe_grouping, x = "Functional.Group", y = "Eigenvalue", facet.by = "supp") +
#   stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)
# 

# devEMF::emf(file = "Anova_emmeans_FGs_violinplot.emf", emfPlus = F)
bxp_emm <- ggviolin(bxp_df, x = "Functional.Group", y = "eigenvalue", 
                    facet.by = "dimension",
                    panel.labs = list(dimension = c("Dimension 1 (43.8%)", 
                                                    "Dimension 2 (33.7%)", 
                                                    "Dimension 3 (22.4%)")),
                    color = "Functional.Group",
                    # fill = "gray", 
                    # palette = c('#1972A4', '#965F8A', '#FF7070', '#4AC6B7'),
                    add = "jitter", bxp.errorbar=T,
                    ylab = "Eigenvalue", xlab = "Functional Groups")
ggpar(bxp_emm, legend = "right") +
  stat_pvalue_manual(bxp_data, hide.ns = T, #Use significant p.values
                     y.position = 0.75,
                     step.increase = 0.1,
                     label = "p.adj.signif",
                     color = "black",
                     size = 5) 

### Coefficient of Wetness is ordinal... 
# spe_grouping$CoW_o <- ordered(as.numeric(spe_grouping$Coefficient.of.Wetness))
scatter_df <- spe_grouping %>% pivot_longer(cols = c(V1:V3), #Create new columns in long form for surveyor ids
                              names_to = c("dimension"),
                              values_to = "eigenvalues")

#### Scatter Plot ####
sct_emm <- ggplot(scatter_df, aes(x = Coefficient.of.Wetness, y = eigenvalues, 
                                  shape = dimension, color = dimension),
                    panel.labs = list(dimension = c("Dimension 1 (43.8%)",
                                                    "Dimension 2 (33.7%)",
                                                    "Dimension 3 (22.4%)")),
                    # color = "CoW",
                    # fill = "gray",
                    # palette = c('#1972A4', '#965F8A', '#FF7070', '#4AC6B7'),
                  add = "jitter",
                  ylab = "Eigenvalue", xlab = "Coefficient of Wetness") +
  geom_point(size=2) +
  geom_smooth(method="lm", se=TRUE, level=0.95, aes(linetype = dimension)) + 
  labs(x = "Coefficient of Wetness",
       y = "Eigenvalues",
       )

sct_emm + theme_classic() + 
  facet_grid( ~dimension) +
  scale_color_brewer(palette="Dark2") + 
  scale_shape_manual(values=c(3, 16, 17))


# dev.off()




plot(V1~Functional.Group,spe_grouping)
plot(V2~Functional.Group,spe_grouping)
plot(V3~Functional.Group,spe_grouping)

spe.h.cmd$species

x <- data.frame(V1 = spe_grouping$V1, 
                V2 = spe_grouping$V2, 
                V3 = spe_grouping$V3)

rownames(x) <- rownames(spe.h.cmd$species)


# devEMF::emf(file = "LinePlots_speciesV2.emf", emfPlus = F)
par(mfrow = c(3,1))
plot(0, xlim = c(-1, 1), axes=FALSE,
     type = "n", xlab = "Dimension 1 (43.8%)", ylab = "")
axis(1, at = x1, labels = names(x1), par(las = 2))

plot(0, xlim = c(-1, 1), axes=FALSE,
     type = "n", xlab = "Dimension 2 (33.7%)", ylab = "")
axis(1, at = x2, labels = names(x2))

plot(0, xlim = c(-1, 1), axes=FALSE,
     type = "n", xlab = "Dimension 3 (22.4%)", ylab = "")
axis(1, at = x3, labels = names(x3))
# dev.off()

library(ggpubr)

# lp1 <- ggplot(data = x, aes(x = V1, y = 0)) + 
#   geom_blank() + theme_classic() + 
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_x_continuous(breaks = x$V1, labels = rownames(x), 
#                      guide = guide_axis(n.dodge=1), #Text dodge function (to play around with)
#                      name = "Dimension 1 (43.8%)") +
#   scale_y_continuous(name = "")
# lp2 <- ggplot(data = x, aes(x = V2, 0)) + 
#   geom_blank() + theme_classic() + 
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_x_continuous(breaks = x$V2, labels = rownames(x), 
#                      name = "Dimension 2 (33.7%)")+
#   scale_y_continuous(name = "")
# lp3 <- ggplot(data = x, aes(x = V3, 0)) + 
#   geom_blank() + theme_classic() + 
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_x_continuous(breaks = x$V3, labels = rownames(x), 
#                      name = "Dimension 3 (22.4%)")+
#   scale_y_continuous(name = "")
# 
# line_graphs <- ggarrange(lp1, lp2, lp3, 
#                          ncol = 1, nrow = 3)
# fc_graphs <- ggarrange(fc_12, fc_13, fc_23,
#                        nrow = 2, ncol = 2)

devEMF::emf(file = "LinePlots_species.emf", emfPlus = F)

dev.off()

#3D cluster plots?
library(tidyverse)
library(plotly)
# library(htmlwidgets)
colors_p <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951', 'black')
plot_ly(as.data.frame(spe.h.cmd$species), x = ~ V1, 
        y = ~ V2, 
        z = ~ V3, 
        type = "scatter3d", 
        mode = "markers", 
        color = as.factor(final$clustering),
        colors = colors_p
) %>%
  add_markers() %>% 
  layout(title = "",
         scene = list(xaxis = list(title = "Dimension 1 (43.8%)"),
                      yaxis = list(title = "Dimension 2 (33.7%)"),
                      zaxis = list(title = "Dimension 3 (22.4%)")
         )
  )


plot(stem(x = spe.h.cmd$species[,1]))


library(scatterplot3d)
scatterplot3d(spe.h.cmd$species[,1:3], pch=final$clustering,
              color = c('purple', 'black', 'blue', 'darkgreen', 'grey25', 'brown')[final$clustering],
              angle = 20, 
              tick.marks = F, 
              xlab= "Dim 1", ylab="Dim 2", zlab="Dim 3")

ggplot(data=spe_grouping,aes(x=V1,y=V3,col=Functional.Group)) + geom_point(size=5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw()


###############
library(ggordiplots)
ordiplot(spe.h.cmd, choices = c(1,2), type = "text",
         cex = 0.6, cex.axis = .6, ann = F, tck=-0.01) 
plot(envfit(spe.h.cmd, spe.dt[,-c(1, 257:262)])) ##### Got it!

plt1 <- gg_ordiplot(spe.h.cmd$points, groups = tmptest$Functional.Group, plot = FALSE)

wascores(spe.h.cmd$points, tmptest[,-c(1, 257:262)], expand = T)

ordiplot(spe.h.cmd, choices = c(1,2), type = "text",
         cex = 0.6, cex.axis = .6, ann = F, tck=-0.01) 

ordiplot(spe.h.cmd, display = c('sites'), type = 'p') 
orditorp(spe.h.cmd, display = 'sp')
ordiplot(spe.h.cmd, type = 'point')

#Use the spe.h.cmd$species dataset to combine the functional groups
rownames(spe.h.cmd$species)

library (vegan3d)
ordirgl(spe.h.cmd)

# plot(envfit(spe.h.cmd, tmp4[,-c(1)])) #RRA_ids and rows from points are not matching

###################

ordiplot(spe.h.cmd, choices = c(1,3), type = "text",
       cex = 0.6, cex.axis = .6, ann = F, tck=-0.01) 
### Select only the top 20 species to plot? How to do this
plot(envfit(spe.h.cmd, spe.dt[,-c(1, 257:2262)])) ##### Got it! 



library(ggordiplots)
gg_envfit(spe.h.cmd)



spe.dt$cmd1 <- spe.h.cmd$points[,1]  
spe.dt$cmd2 <- spe.h.cmd$points[,2]
spe.dt$cmd3 <- spe.h.cmd$points[,3]

#### Surveyor NMDS Supplementary ####
# devEMF::emf(file='NMDS_surveyors.emf',emfPlus = F)
par(mfrow=c(2,2))
#A - Stress plot across multiple dimensions
barplot(stress, ylab = "stress", 
        main = "A") 
abline(h = 0.01, col = "darkred")
#B #Shepard plot
spear <- round(cor(vegdist(surv, method = "jaccard"), 
                   dist(surv.nmds$points), 
                   method = "spearman"),3) #Spearman test
stressplot(surv.nmds, main="B"
) 
text(x = 0.2, y =  7.20, #Hardcode the position of the spearman test text
     paste0("Spearman correlation = ", spear))
#C #paste("Goodness of Fit - NMDS/Jaccard - Stress =", round(surv.nmds$stress,3))
eig_1 <- round(surv.cmd$eig*100/sum(surv.cmd$eig),1)[1]
eig_2 <- round(surv.cmd$eig*100/sum(surv.cmd$eig),1)[2]
plot(surv.cmd$points, type = "n", 
     main = "C",
     xlab = paste0("MDS-1 (", eig_1, "%)"),
     ylab = paste0("MDS-2 (", eig_2, "%)")
)
points(surv.cmd$points, cex = 2, pch = 21, col = "black", bg = 'darkred') #Poorly fit sites have large bubbles
orditorp(surv.cmd, display = "species", 
         air = 1, cex = 1) #Surveyor labels are given space for readability 
#D
gof <- goodness(surv.nmds)
plot(surv.nmds$points, type = "n", 
     main = "D", #paste("Goodness of Fit - NMDS/Jaccard - Stress =", round(surv.nmds$stress,3))
     # xlab = paste0("MDS-1 (", eig_1, "%)"),
     # ylab = paste0("MDS-2 (", eig_2, "%)")
)
points(surv.nmds$points, pch = 21, col = "black", bg = 'steelblue',
       cex = gof*2/mean(gof) #Poorly fit sites have large bubbles
)
orditorp(surv.nmds, display = "species", 
         air = 1, cex = 1 #Surveyor labels are given space for readability
)
# dev.off()




#S2
library(sjPlot)

#S3
sjPlot::tab_model(min.HerbInt, min.BExt,
                  digits = 6,
                  #    pred.labels =c("(Intercept) - 10% Survival", "(Intercept) - 25% Survival", 
                  # "(Intercept) - 50% Survival", "(Intercept) - 75% Survival",
                  # "Simpson diversity", "Community composition Axis #1", 
                  # "Community composition Axis #2", "Community composition Axis #3",
                  # "Urbanization (PC1)", "Impervious land uses (PC2)", "Distance to city centre",
                  # "Planting month",
                  # "Human interferrence richness", "Herbivore richness", "Invasive species richness", 
                  # "Assessment year", "Surveyor bias"),
                  dv.labels= c("Herbivory intensity", "Herbivory extent"), 
                  show.intercept = F,
                  show.df = F #Inf? why
                  
)
anova(cm.survival, cn.survival) #fixed effects 
anova(cm.health.min, cn.health)

#S4
sjPlot::tab_model(min.InvSpp, min.InvExt,
                  digits = 6,
                  #    pred.labels =c("(Intercept) - 10% Survival", "(Intercept) - 25% Survival", 
                  # "(Intercept) - 50% Survival", "(Intercept) - 75% Survival",
                  # "Simpson diversity", "Community composition Axis #1", 
                  # "Community composition Axis #2", "Community composition Axis #3",
                  # "Urbanization (PC1)", "Impervious land uses (PC2)", "Distance to city centre",
                  # "Planting month",
                  # "Human interferrence richness", "Herbivore richness", "Invasive species richness", 
                  # "Assessment year", "Surveyor bias"),
                  dv.labels= c("Invasive species intensity", "Invasive species extent"), 
                  show.intercept = F,
                  show.df = F #Inf? why
                  
)

#S5

sjPlot::tab_model(min.HIS, min.HumanExt,
                  digits = 6,
                  #    pred.labels =c("(Intercept) - 10% Survival", "(Intercept) - 25% Survival", 
                  # "(Intercept) - 50% Survival", "(Intercept) - 75% Survival",
                  # "Simpson diversity", "Community composition Axis #1", 
                  # "Community composition Axis #2", "Community composition Axis #3",
                  # "Urbanization (PC1)", "Impervious land uses (PC2)", "Distance to city centre",
                  # "Planting month",
                  # "Human interferrence richness", "Herbivore richness", "Invasive species richness", 
                  # "Assessment year", "Surveyor bias"),
                  dv.labels= c("Human disturbance intensity", "Human disturbance extent"), 
                  show.intercept = F,
                  show.df = F #Inf? why
                  
)



# fixed effects model
t1.survival <- clm(Survival_score ~  hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                      pc1 + pc2 + 
                     # dist_TO + #Landscape Factors
                      planting_month + #Dispersal/phenology
                      NumHuman_Interf + 
                     # browsing_rich + 
                     Obs_invSpecies + #Disturbances
                      # Assessment.Year. + 
                     planting_type, nominal = ~ pcoa_surv1,  #Survey
                    Hess = T,
                    data = clmm_surv[!is.na(clmm_surv$Survival_score),])
# mixed model
cm.survival <- clmm(Survival_score ~  hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                      pc1 + pc2 + dist_TO + #Landscape Factors
                      planting_month + #Dispersal/phenology 
                      NumHuman_Interf + browsing_rich + Obs_invSpecies + #Disturbances
                      Assessment.Year. + planting_type + pcoa_surv1 +#Survey
                      (1|RRA_id), #Survey as random variable.
                    Hess = T,
                    data = clmm_surv[!is.na(clmm_surv$Survival_score),])



#Landscape 
ggcorrplot(RRA_LU_cor) #Landuse correlation plot of reclassed landuses
summary(RRA_LU_pca) #Summary table of PCA and cumulative contributions


#Herbivory
library(ggridges)
clmm_surv %>% 
  ggplot(aes(x = browse_deterrents,
             y = Browsing.of.planted.vegitation.,
             fill = num_browsDeterrents)) +
  ggridges::geom_density_ridges(bandwidth = 0.5)

test_browse <- clm(Browsing.of.planted.vegitation. ~ browse_deterrents, 
                   Hess = T,
                   data = clmm_surv[!is.na(clmm_surv$Browsing.of.planted.vegitation.),])
summary(test_browse)
anova(test_browse, type = "III")



################################################################################

### Herbivory ###
#Herbivory = 3 heavy browsing < 0 no browsing
lm.HerbInt <- lm(as.numeric(fct_rev(Browsing.of.planted.vegitation.)) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                   pc1 + pc2 + dist_TO + #Landscape Factors
                   planting_month + #Restoration  
                   Assessment.Year. + pcoa_surv1 +
                   browse_deterrents + NumHuman_Interf + Obs_invSpecies,
                 data = clmm_data[!is.na(clmm_data$Browsing.of.planted.vegitation.),])
summary(lm.HerbInt)
vif(lm.HerbInt)

# cn.HerbInt <- clmm(fct_rev(Browsing.of.planted.vegitation.) ~ (1|RRA_id), 
#                data = clmm_surv[!is.na(clmm_surv$Browsing.of.planted.vegitation.),])

fix.HerbInt <- clm(fct_rev(Browsing.of.planted.vegitation.) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                     pc1 + pc2 + dist_TO + #Landscape Factors
                     planting_month + #Restoration  
                     Assessment.Year. + pcoa_surv1 +
                     browsing_rich + 
                     num_browsDeterrents,
                   Hess = T, 
                   data = clmm_data[!is.na(clmm_data$Browsing.of.planted.vegitation.),])

# cm.HerbInt <- clmm(Browsing.of.planted.vegitation. ~ hill_shan + pcoa1 + pcoa2 + pcoa3 + #Planting Diversity and Composition
#                      pc1 +
#                      pc2 +
#                      dist_TO + #Landscape Factors
#                      # planting_month + #Restoration
#                      Assessment.Year. +
#                      pcoa_surv1 +
#                      browse_deterrents +
#                      (1|RRA_id),
#                   Hess = T,
#                   data = clmm_surv[!is.na(clmm_surv$Browsing.of.planted.vegitation.),])

# summary(cn.HerbInt)
summary(fix.HerbInt)
# summary(cm.HerbInt)

anova(fix.HerbInt, type = "II")
# anova(cm.HerbInt, cn.HerbInt)
# anova(cm.HerbInt, fix.HerbInt) 

stepf_HerbInt <- step(fix.HerbInt, direction = "both")
# Minimal model
min.HerbInt <- clm(fct_rev(Browsing.of.planted.vegitation.) ~  cmd1 + cmd3 + pc1 + 
                     pc2 + pcoa_surv1 + browsing_rich + browse_deterrents,
                   Hess = T,
                   data = clmm_data[!is.na(clmm_data$Browsing.of.planted.vegitation.),])
summary(min.HerbInt)

## Extent of Browsing ##
fix.BExt <- clm(fct_rev(Extent.of.browsing.) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                  pc1 + pc2 + dist_TO + #Landscape Factors
                  planting_month + #Restoration
                  Assessment.Year. +
                  pcoa_surv1 + browsing_rich +
                  browse_deterrents, 
                Hess = T, 
                data = clmm_data[!is.na(clmm_data$Extent.of.browsing.),])


summary(fix.BExt)
anova(fix.BExt, type = "II") #***
# Minimal model
stepf_BExt <- step(fix.BExt, direction = "both")
min.BExt <- clm(fct_rev(Extent.of.browsing.) ~  cmd1 + cmd2 + pc2 + pcoa_surv1 + 
                  browsing_rich + browse_deterrents,
                Hess = T,
                data = clmm_data[!is.na(clmm_data$Extent.of.browsing.),])
summary(min.BExt)


### Invasive Species ###
lm.InvSpp <- lm(as.numeric(fct_rev(Invasive.species.)) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                  pc1 + pc2 + dist_TO + #Landscape Factors
                  planting_month + #Restoration  
                  Assessment.Year. + pcoa_surv1 +
                  Obs_invSpecies + Year.last.treated.,
                data = clmm_data[!is.na(clmm_data$Invasive.species.),])
summary(lm.InvSpp)
vif(lm.InvSpp)


#### Invasive Species 3 bad < 0 good

fix.InvSpp <- clm(fct_rev(Invasive.species.) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                    pc1 + pc2 + dist_TO + #Landscape Factors
                    Assessment.Year. +
                    pcoa_surv1 + Obs_invSpecies +
                    Year.last.treated., 
                  Hess = T, 
                  data = clmm_data[!is.na(clmm_data$Invasive.species.),])

summary(fix.InvSpp)
anova(fix.InvSpp, type = "II")

# Minimal model
stepf_InvSpp <- step(fix.InvSpp, direction = "both")
min.InvSpp <- clm(fct_rev(Invasive.species.) ~  cmd1 + pc1 + pc2 + pcoa_surv1,
                  Hess = T,
                  data = clmm_data[!is.na(clmm_data$Invasive.species.),])
summary(min.InvSpp)


### Extent of Invasive species spread ##

fix.InvExt <- clm(fct_rev(Extent.of.invasive.species.) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                    pc1 + pc2 + dist_TO + #Landscape Factors
                    Assessment.Year. +
                    pcoa_surv1 + Obs_invSpecies +
                    Year.last.treated., 
                  Hess = T, 
                  data = clmm_data[!is.na(clmm_data$Invasive.species.),])
summary(fix.InvExt)
# Minimal model
stepf_InvExt <- step(fix.InvExt, direction = "both")
min.InvExt <- clm(fct_rev(Extent.of.browsing.) ~  pc1 + pc2 + pcoa_surv1 + 
                    Obs_invSpecies + Year.last.treated.,
                  Hess = T,
                  data = clmm_data[!is.na(clmm_data$Extent.of.browsing.),])
summary(min.InvExt)

#### Human Disturbance ####
fix.HIS <- clm(fct_rev(Human.interference.with.site.) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                 pc1 + pc2 + dist_TO + #Landscape Factors
                 Assessment.Year. +
                 pcoa_surv1 + 
                 NumHuman_Interf + NumHuman_deterrents,
               Hess = T, 
               data = clmm_data[!is.na(clmm_data$Human.interference.with.site.),]) 

summary(fix.HIS)
# Minimal model
stepf_HIS <- step(fix.HIS, direction = "both")
min.HIS <- clm(fct_rev(Human.interference.with.site.) ~  cmd1 + cmd3 + pc2 + 
                 pcoa_surv1 + NumHuman_Interf + NumHuman_deterrents,
               Hess = T,
               data = clmm_data[!is.na(clmm_data$Human.interference.with.site.),])
summary(min.HIS)

### Extent of human disturbance ###
fix.HumanExt <- clm(fct_rev(Extent.of.human.interference.) ~ hill_shan + cmd1 + cmd2 + cmd3 + #Planting Diversity and Composition
                      pc1 + pc2 + dist_TO + #Landscape Factors
                      Assessment.Year. +
                      pcoa_surv1 + 
                      NumHuman_Interf + NumHuman_deterrents,
                    Hess = T, 
                    data = clmm_data[!is.na(clmm_data$Human.interference.with.site.),]) 

summary(fix.HumanExt)
# Minimal model
stepf_HumanExt <- step(fix.HumanExt, direction = "both")
min.HumanExt <- clm(fct_rev(Human.interference.with.site.) ~  cmd1 + cmd2 + cmd3 + 
                      pc2 + Assessment.Year. + NumHuman_Interf,
                    Hess = T,
                    data = clmm_data[!is.na(clmm_data$Human.interference.with.site.),])
summary(min.HumanExt)



## Disturbance Intensity ##
# Herbivory #
#Herbivory intensity
herbivory_plot <- plot_model(min.HerbInt,
                             terms = c("pcoa1", "pcoa3",
                                       "pc1","pc2",
                                       "pcoa_surv1",
                                       "browsing_rich",
                                       "browse_deterrentsbeaver guards guards rodent",
                                       "browse_deterrentsdeer fence",
                                       "browse_deterrentsdeer fence guards rodent",
                                       "browse_deterrentsguards rodent"),
                             auto.label = T,
                             vline.color = "blue",
                             show.intercept = F,
                             show.values=TRUE,
                             value.offset = .3,
                             show.p=TRUE,
                             line.size = 1,
                             dot.size = 3,
                             title="",
                             show.data = T)
devEMF::emf(file='Figure_HerbInt.emf',emfPlus = F)
herbivory_plot
dev.off()

#Herbivory extent
# Herbivory #
herbExt_plot <- plot_model(min.BExt,
                           terms = c("pcoa1", "pcoa3",
                                     "pc2","pcoa_surv1",
                                     "browsing_rich","browse_deterrentsbeaver guards guards rodent",
                                     "browse_deterrentsdeer fence",
                                     "browse_deterrentsdeer fence guards rodent",
                                     "browse_deterrentsguards rodent"),
                           auto.label = T,
                           vline.color = "blue",
                           show.intercept = F,
                           show.values=TRUE,
                           value.offset = .3,
                           show.p=TRUE,
                           line.size = 1,
                           dot.size = 3,
                           title="",
                           show.data = T)
devEMF::emf(file='Figure_HerbExt.emf',emfPlus = F)
herbExt_plot
dev.off()


# Human Disturbance #
#Human intensity
human_plot <- plot_model(min.HIS,
                         terms = c("pcoa1", "pcoa3",
                                   "pc2",
                                   "pcoa_surv1",
                                   "NumHuman_Interf",
                                   "NumHuman_deterrents"),
                         auto.label = T,
                         vline.color = "blue",
                         show.intercept = F,
                         show.values=TRUE,
                         value.offset = .3,
                         show.p=TRUE,
                         line.size = 1,
                         dot.size = 3,
                         title="",
                         show.data = T)
devEMF::emf(file='Figure_HumanInt.emf',emfPlus = F)
human_plot
dev.off()

#Human intensity
humanExt_plot <- plot_model(min.HumanExt,
                            terms = c("pcoa1", "pcoa2", "pcoa3",
                                      "pc2",
                                      "Assessment.Year.",
                                      "NumHuman_Interf"),
                            auto.label = T,
                            vline.color = "blue",
                            show.intercept = F,
                            show.values=TRUE,
                            value.offset = .3,
                            show.p=TRUE,
                            line.size = 1,
                            dot.size = 3,
                            title="",
                            show.data = T)
devEMF::emf(file='Figure_HumanExt.emf',emfPlus = F)
humanExt_plot
dev.off()



# Invasive Disturbance #
#Invasive intensity
invasive_plot <- plot_model(min.InvSpp,
                            terms = c("pcoa1",
                                      "pc1", "pc2",
                                      "pcoa_surv1"),
                            auto.label = T,
                            vline.color = "blue",
                            show.intercept = F,
                            show.values=TRUE,
                            value.offset = .3,
                            show.p=TRUE,
                            line.size = 1,
                            dot.size = 3,
                            title="",
                            show.data = T)
devEMF::emf(file='Figure_InvasiveInt.emf',emfPlus = F)
invasive_plot
dev.off()

#Invasive Extent
invasiveExt_plot <- plot_model(min.InvExt,
                               terms = c("pc1", "pc2",
                                         "pcoa_surv1", 
                                         "Obs_invSpecies",
                                         "Year.last.treated."),
                               auto.label = T,
                               vline.color = "blue",
                               show.intercept = F,
                               show.values=TRUE,
                               value.offset = .3,
                               show.p=TRUE,
                               line.size = 1,
                               dot.size = 3,
                               title="",
                               show.data = T)
devEMF::emf(file='Figure_InvasiveExt.emf',emfPlus = F)
invasiveExt_plot
dev.off()
