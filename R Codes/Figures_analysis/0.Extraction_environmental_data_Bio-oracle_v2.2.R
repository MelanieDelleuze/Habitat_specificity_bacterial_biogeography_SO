#### Import environmental data from BIO-Oracle v.2.2 ####
# Assis et al. 2017
# https://bio-oracle.org/index.php

#Set working directory
setwd("~/GitHub/Habitat_specificity_bacterial_biogeography_SO")

# Only once and then save the metadata with the extracted environmental variables (phyloseq@samdata)
library(sdmpredictors)
### List of layers available
layers.bio2 = list_layers(datasets="Bio-ORACLE" )
View(layers.bio2)
write.table(layers.bio2, file = "bio_oracle.layers.txt", sep = "\t",row.names = TRUE, col.names = NA)

### Dataframe containing information of the sites and samples
#metadata phyloseq with a column with latitude and another with longitude for each sample
my_data = read.table("metadata.habitat.specificity.SO.txt",header=TRUE) 

# Download the environmental layers of interest
options(timeout = 1000) #increasing timeout
layers.envi =load_layers(layercodes = c("BO22_chlomean_bdmean",
"BO22_curvelmean_bdmean","BO22_dissoxmean_bdmean","BO22_ironmean_bdmean",
"BO22_phosphatemean_bdmean","BO22_lightbotmean_bdmean","BO22_ph","BO22_nitratemean_bdmean",
"BO22_tempmean_bdmean","BO22_carbonphytomean_bdmean","BO22_ppmean_bdmean",
"BO22_salinitymean_bdmean","BO22_silicatemean_bdmean"),equalarea = FALSE, rasterstack = TRUE,datadir ="~/GitHub/Habitat_specificity_bacterial_biogeography_SO/Bio_Oracle")

# Extracting the environmental values for each of the sample's site
my_data_environment = data.frame(Name=my_data$Group ,raster::extract(layers.envi,my_data[,14:15])) #lat and long columns

# Scale (center reduce) environmental data
my_data_environment_cr = scale(my_data_environment[,2:14],center = TRUE,scale = TRUE)
# Add the environmental data to the initial metadata
my_data_envi = cbind(my_data[,1:15],my_data_environment_cr)

# Saving the metadata with environmental variables, this is the metadata used for the phyloseq object
write.table(my_data_envi, file = "metadata.habitat.specificity.txt", sep = "\t",row.names = TRUE, col.names = NA) 




