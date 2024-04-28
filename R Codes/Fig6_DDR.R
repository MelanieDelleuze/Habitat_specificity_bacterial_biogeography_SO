####------------------------------------------------##
#### Fig.6: Geographical and Environmental DDR    ####
##--------------------------------------------------##

#Load packages
library(phyloseq)
library(microViz)
library(geosphere)   
library(reshape) 
library(vegan)
library(ggplot2)
library(geosphere)        
library(usdm)
library(dplyr)
library(ggpubr)
library(cowplot)
library(tidyverse)

#Set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#Load data
load("physeq_habitat_specificity.RData") 
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
content = ps_filter(physeq_habitat_specificity,Type=="Gut content")
tissue = ps_filter(physeq_habitat_specificity,Type=="Gut tissue")
# Transform phyloseqs in relative abundance
sediment.rel = microbiome::transform(sediment,"compositional")
content.rel = microbiome::transform(content,"compositional")
tissue.rel = microbiome::transform(tissue,"compositional")
#Bray-Curtis dissimilarity matrix
dist.sed = phyloseq::distance(sediment.rel, method="bray")
dist.cont = phyloseq::distance(content.rel, method="bray")
dist.tiss = phyloseq::distance(tissue.rel, method="bray")
#Linearise Bray-Curtis matrix
bc.sed=melt(as.matrix(dist.sed))
bc.cont=melt(as.matrix(dist.cont))
bc.tiss=melt(as.matrix(dist.tiss))

#### Geographic distances ####
# Create a dataframe with the longitude, latitude and samples for each habitat
# Sediment
long.sed = sediment.rel@sam_data$Longitude
lat.sed = sediment.rel@sam_data$Latitude
sample.sed = sediment.rel@sam_data$Group
lats.sed = data.frame(long.sed,lat.sed,sample.sed)
rownames(lats.sed) = lats.sed$sample.sed
lats.sed2= lats.sed[,-3]
#Gut content
long.cont = content.rel@sam_data$Longitude
lat.cont = content.rel@sam_data$Latitude
sample.cont = content.rel@sam_data$Group
lats.cont = data.frame(long.cont,lat.cont,sample.cont)
rownames(lats.cont) = lats.cont$sample.cont
lats.cont2= lats.cont[,-3]
#Gut tissue
long.tiss = tissue.rel@sam_data$Longitude
lat.tiss = tissue.rel@sam_data$Latitude
sample.tiss = tissue.rel@sam_data$Group
lats.tiss = data.frame(long.tiss,lat.tiss,sample.tiss)
rownames(lats.tiss) = lats.tiss$sample.tiss
lats.tiss2= lats.tiss[,-3]

# Calculate the distances between sites in meters in a ellipsoid map WGS-84
#(distGeo function does not provide a symmetric matrix as distm function) -> loop to do it

# Initialize an empty matrix to store distances
data.coord =  lats.tiss2   # lats.sed2 / lats.cont2 /lats.tiss2 #change by your table with longitude and latitude
n <- nrow(data.coord)
distance_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(distance_matrix) <- rownames(distance_matrix) <- rownames(data.coord)
# Calculate distances between all pairs of samples
for (i in 1:n) {
  for (j in 1:n) {
    # Extract latitude and longitude for samples i and j
    coords_i <- as.numeric(unlist(data.coord[i, ]))
    coords_j <- as.numeric(unlist(data.coord[j, ]))
    # Calculate distance using distGeo
    distance_matrix[i, j] <- distGeo(coords_i, coords_j)
    distance_matrix[j, i] <- distance_matrix[i, j] # Make the matrix symmetric
  }
}
#Store the matrix in other object for each one of the habitats
#dist_geo_sed=distance_matrix
#dist_geo_cont=distance_matrix
dist_geo_tiss=distance_matrix

# Linearise distance matrix for later (construct data frame with all distances)
geo.sed=melt(as.matrix(dist_geo_sed))
geo.cont=melt(as.matrix(dist_geo_cont))
geo.tiss=melt(as.matrix(dist_geo_tiss))

# Transform to km and linearise
dist_geo_sed_km = dist_geo_sed/1000
geo.sed.km = melt(as.matrix(dist_geo_sed_km))

dist_geo_cont_km = dist_geo_cont/1000
geo.cont.km = melt(as.matrix(dist_geo_cont_km))

dist_geo_tiss_km = dist_geo_tiss/1000
geo.tiss.km = melt(as.matrix(dist_geo_tiss_km))

#### Environmental distance ####
# Sediment
data_envi_sed0 = data.frame(sediment.rel@sam_data)
data_envi_sed = data_envi_sed0[,-(1:15)]
#13 environmental variables

#Keeping the same environmental variables as the ones represented in PCA
data_envi_sed_1 = data_envi_sed[c("BO22_curvelmean_bdmean", "BO22_dissoxmean_bdmean",
                                "BO22_ironmean_bdmean" , "BO22_lightbotmean_bdmean"   ,
                                "BO22_ph" , "BO22_nitratemean_bdmean" ,    "BO22_tempmean_bdmean" ,      
                                "BO22_salinitymean_bdmean"   ,
                                "BO22_silicatemean_bdmean" )] #9 environmental variables

# Creating a environmental distance matrix between samples with data_envi_sed_1 (without high correlated ones)
dist.env.sed = dist(data_envi_sed_1, method = "euclidean")
env.sed = melt(as.matrix(dist.env.sed)) # Linearise the matrix -> data frame for later


# Gut content
data_envi_cont0 = data.frame(content.rel@sam_data)
data_envi_cont = data_envi_cont0[,-(1:15)]
#13 environmental variables

#Keeping the same environmental variables as the ones represented in PCA
data_envi_cont_1 = data_envi_cont[c("BO22_curvelmean_bdmean", "BO22_dissoxmean_bdmean",
                                  "BO22_ironmean_bdmean" , "BO22_lightbotmean_bdmean"   ,
                                  "BO22_ph" , "BO22_nitratemean_bdmean" ,    "BO22_tempmean_bdmean" ,      
                                  "BO22_salinitymean_bdmean"   ,
                                  "BO22_silicatemean_bdmean" )] #9 environmental variables

# Creating a environmental distance matrix between samples with data_envi_cont_1 (without high correlated ones)
dist.env.cont = dist(data_envi_cont_1, method = "euclidean")
env.cont = melt(as.matrix(dist.env.cont)) # Linearise the matrix -> data frame for later

# Gut tissue
data_envi_tiss0 = data.frame(tissue.rel@sam_data)
data_envi_tiss = data_envi_tiss0[,-(1:15)]
#13 environmental variables

#Keeping the same environmental variables as the ones represented in PCA
data_envi_tiss_1 = data_envi_tiss[c("BO22_curvelmean_bdmean", "BO22_dissoxmean_bdmean",
                                    "BO22_ironmean_bdmean" , "BO22_lightbotmean_bdmean"   ,
                                    "BO22_ph" , "BO22_nitratemean_bdmean" ,    "BO22_tempmean_bdmean" ,      
                                    "BO22_salinitymean_bdmean"   ,
                                    "BO22_silicatemean_bdmean" )] #9 environmental variables

# Creating a environmental distance matrix between samples with data_envi_tiss_1 (without high correlated ones)
dist.env.tiss = dist(data_envi_tiss_1, method = "euclidean")
env.tiss = melt(as.matrix(dist.env.tiss)) # Linearise the matrix -> data frame for later

#### Creating data frames with all distances for each habitat ####
#Sediment
distances.sed = data.frame(bc.sed,geo.sed$value,geo.sed.km$value,env.sed$value)
colnames(distances.sed) = c("Sample_1","Sample_2","BC_Dissimilarity","Distance_m","Distance_km","Distance_envi")
#Gut content
distances.cont = data.frame(bc.cont,geo.cont$value,geo.cont.km$value,env.cont$value)
colnames(distances.cont) = c("Sample_1","Sample_2","BC_Dissimilarity","Distance_m","Distance_km","Distance_envi")
#Gut tissue
distances.tiss = data.frame(bc.tiss,geo.tiss$value,geo.tiss.km$value,env.tiss$value)
colnames(distances.tiss) = c("Sample_1","Sample_2","BC_Dissimilarity","Distance_m","Distance_km","Distance_envi")

#### DDR ####
##### Sediment #####
# Linear regression model
# Adding + 1 for log transformation and log transformation
distances.sed$dist.km_1= distances.sed$Distance_km + 1
distances.sed$logDist=log10(distances.sed$dist.km_1)

distances.sed$dist.env_1 = distances.sed$Distance_envi + 1
distances.sed$logEnv = log10(distances.sed$dist.env_1)

#Removing same comparisons of samples
distances.sed$Sample_1 = as.character(distances.sed$Sample_1)
distances.sed$Sample_2 = as.character(distances.sed$Sample_2)
distances.sed2 = distances.sed %>% dplyr::filter(Sample_1 < Sample_2)


#fit simple linear regression model
fit <- lm(BC_Dissimilarity ~ logDist, data=distances.sed2) #change logDist by logEnv to have the other model
#view model summary
summary(fit)

# Plot DDR Geo
mantel.plot.sed = ggplot(distances.sed2, aes(y = BC_Dissimilarity, x=logDist)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = Distance_envi)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  stat_regline_equation(label.y = 1.09, aes(label = ..eq.label..),size=13) +
  stat_cor(label.y = 1,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=13)+
  labs(x = 'Geographic distance (km,' ~ log[10]~ ')')+
  labs(y = "Bray-Curtis Dissimilarity", fill = "Geographic distance (km)") + 
  ggtitle("Sediment")+
  scale_x_continuous(limits= c(0,4),breaks = c(0,2,4))+
  scale_y_continuous(limits= c(0,1.1),breaks = c(0,0.5,1))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=38),
        legend.title = element_text(size=38),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))+
  scale_fill_binned(high = "#b90000", low = "#ffdc5e")

mantel.plot.sed

#Plot DDR Envi
mantel.plot.sed.env = ggplot(distances.sed2, aes(y = BC_Dissimilarity, x=logEnv)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = Distance_km)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  stat_regline_equation(label.y = 1.09, aes(label = ..eq.label..),size=13) +
  stat_cor(label.y = 1.02,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=13)+
  labs(x = 'Environmental distance (' ~ log[10]~ ')')+
  labs(y = "Bray-Curtis Dissimilarity", fill = "Geographic distance (km)") + 
  ggtitle("")+
  scale_x_continuous(limits= c(0,1),breaks = c(0,0.5,1))+
  scale_y_continuous(limits= c(0,1.1),breaks = c(0,0.5,1))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=38),
        legend.title = element_text(size=38),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))+
  scale_fill_binned(high = "#b90000", low = "#ffdc5e")

mantel.plot.sed.env

##### Gut content #####
# Linear regression model
# Adding + 1 for log transformation and log transformation
distances.cont$dist.km_1= distances.cont$Distance_km + 1
distances.cont$logDist=log10(distances.cont$dist.km_1)

distances.cont$dist.env_1 = distances.cont$Distance_envi + 1
distances.cont$logEnv = log10(distances.cont$dist.env_1)

#Removing same comparisons of samples
distances.cont$Sample_1 = as.character(distances.cont$Sample_1)
distances.cont$Sample_2 = as.character(distances.cont$Sample_2)
distances.cont2 = distances.cont %>% dplyr::filter(Sample_1 < Sample_2)

#fit simple linear regression model
fit <- lm(BC_Dissimilarity ~ logDist, data=distances.cont) #change logDist by logEnv to have the other model
#view model summary
summary(fit)

# Plot DDR geo
mantel.plot.cont = ggplot(distances.cont2, aes(y = BC_Dissimilarity, x=logDist)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = Distance_envi)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  stat_regline_equation(label.y = 1.09, aes(label = ..eq.label..),size=13) +
  stat_cor(label.y = 1,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=13)+
  labs(x = 'Geographic distance (km,' ~ log[10]~ ')')+
  labs(y = "Bray-Curtis Dissimilarity", fill = "Geographic distance (km)") +   ggtitle("Gut content")+
  scale_x_continuous(limits= c(0,4),breaks = c(0,2,4))+
  scale_y_continuous(limits= c(0,1.1),breaks = c(0,0.5,1))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=38),
        legend.title = element_text(size=38),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))+
  scale_fill_binned(high = "#b90000", low = "#ffdc5e")
mantel.plot.cont

#Plot DDR Envi
mantel.plot.cont.env = ggplot(distances.cont2, aes(y = BC_Dissimilarity, x=logEnv)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = Distance_km)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  stat_regline_equation(label.y = 1.09, aes(label = ..eq.label..),size=13) +
  stat_cor(label.y = 1.02,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=13)+
  labs(x = 'Environmental distance (' ~ log[10]~ ')')+
  labs(y = "Bray-Curtis Dissimilarity", fill = "Geographic distance (km)") + 
  ggtitle("")+
  scale_x_continuous(limits= c(0,1),breaks = c(0,0.5,1))+
  scale_y_continuous(limits= c(0,1.1),breaks = c(0,0.5,1))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=38),
        legend.title = element_text(size=38),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))+
  scale_fill_binned(high = "#b90000", low = "#ffdc5e")
mantel.plot.cont.env


##### Gut tissue #####
# Linear regression model
# Adding + 1 for log transformation and log transformation
distances.tiss$dist.km_1= distances.tiss$Distance_km + 1
distances.tiss$logDist=log10(distances.tiss$dist.km_1)

distances.tiss$dist.env_1 = distances.tiss$Distance_envi + 1
distances.tiss$logEnv = log10(distances.tiss$dist.env_1)

#Removing same comparisons of samples
distances.tiss$Sample_1 = as.character(distances.tiss$Sample_1)
distances.tiss$Sample_2 = as.character(distances.tiss$Sample_2)
distances.tiss2 = distances.tiss %>% dplyr::filter(Sample_1 < Sample_2)

#fit simple linear regression model
fit <- lm(BC_Dissimilarity ~ logDist, data=distances.tiss) #change logDist by logEnv to have the other model
#view model summary
summary(fit)


# Plot DDR Geo
mantel.plot.tiss = ggplot(distances.tiss2, aes(y = BC_Dissimilarity, x=logDist)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = Distance_envi)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  stat_regline_equation(label.y = 1.09, aes(label = ..eq.label..),size=13) +
  stat_cor(label.y = 1,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=13)+
  labs(x = 'Geographic distance (km,' ~ log[10]~ ')')+
  labs(y = "Bray-Curtis Dissimilarity", fill = "Environmental distance") +   ggtitle("Gut tissue")+
  scale_x_continuous(limits= c(0,4),breaks = c(0,2,4))+
  scale_y_continuous(limits= c(0,1.1),breaks = c(0,0.5,1))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.position="none",
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(5, 'cm'),
        legend.text = element_text(size=44,  vjust = 0.5),
        legend.title = element_text(size=44, hjust=0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))+
  scale_fill_binned(high = "#b90000", low = "#ffdc5e")

mantel.plot.tiss

#Plot DDR Envi
mantel.plot.tiss.env = ggplot(distances.tiss2, aes(y = BC_Dissimilarity, x=logEnv)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = Distance_km)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  stat_regline_equation(label.y = 1.09, aes(label = ..eq.label..),size=13) +
  stat_cor(label.y = 1.02,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=13)+
  labs(x = 'Environmental distance (' ~ log[10]~ ')')+
  labs(y = "Bray-Curtis Dissimilarity", fill = "Geographic distance (km)") + 
  ggtitle("")+guides(size = guide_legend(title.position="top", title.hjust = 0.5))+
  scale_x_continuous(limits= c(0,1),breaks = c(0,0.5,1))+
  scale_y_continuous(limits= c(0,1.1),breaks = c(0,0.5,1))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.position="none",
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(5, 'cm'),
        legend.text = element_text(size=44, vjust = 0.5),
        legend.title = element_text(size=44, hjust=0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))+
  scale_fill_binned(high = "#b90000", low = "#ffdc5e")

mantel.plot.tiss.env

#Extract legends
leg.ddr.geo <-  ggpubr::get_legend(mantel.plot.tiss)
leg.ddr.geo2 = as_ggplot(leg.ddr.geo)
#save tiff 1000 500

leg.ddr.env <- ggpubr::get_legend(mantel.plot.tiss.env)
leg.ddr.env2 = as_ggplot(leg.ddr.env)
#save tiff 1000 500

#### Bootstrap DDR slopes ####
#Function to subset the phyloseq
sample_ps <- function(ps, FUN = sample, ...){
  ids <- sample_names(ps)
  sampled_ids <- FUN(ids, ...)
  ps <- prune_samples(sampled_ids, ps)
  return(ps)
}

# Bootstrap function to resample 21 samples and calculate DDR slopes
ps=sediment #change phyloseq
iterations=1000 #number of bootstraps
p=21 #number of samples to bootstraps here 75% of 28 = 21

slopes_env <- numeric(iterations) #vector to save the data
slopes_geo <- numeric(iterations)

phylo_samples <- data.frame(Samples=character(), iterations=numeric())

slopes_data_frame <- data.frame(iterations=numeric(), slope_ddr_env=numeric(), slope_ddr_geo=numeric())

set.seed(1234)
for (i in 1:iterations) {
  # Sample 21 samples
  ps_sampled =sample_ps(ps, size=p)
  #Save samples names in a dataframe
  samples_names= ps_sampled@sam_data$Group
  phylo_samples=rbind(phylo_samples, data.frame(Samples=samples_names, iterations=i))
  
  # Calculate distances
  # Bray curtis
  dist = phyloseq::distance(ps_sampled, method="bray")
  bc_dist = melt(as.matrix(dist))
  
  # Geographic data
  sam_data_ps = data.frame(sample_data(ps_sampled))
  sam_data_ps_geo = sam_data_ps[c("Longitude","Latitude")]
  
  # Calculate geographic distances with distGeo function
  data.coord =  sam_data_ps_geo
  n <- nrow(data.coord)
  distance_matrix <- matrix(NA, nrow = n, ncol = n)
  colnames(distance_matrix) <- rownames(distance_matrix) <- rownames(data.coord)
  
  for (x in 1:n) {
    for (j in 1:n) {
      coords_i <- as.numeric(unlist(data.coord[x, ]))
      coords_j <- as.numeric(unlist(data.coord[j, ]))
      distance_matrix[x, j] <- distGeo(coords_i, coords_j)
      distance_matrix[j, x] <- distance_matrix[x, j]
    }
  }
  dist_geo_km = distance_matrix/1000 #transform to km
  geo=melt(as.matrix(dist_geo_km)) #linearise
  
  # Environmental data
  data_envi = data.frame(ps_sampled@sam_data)
  data_envi1= data_envi[,-(1:15)]
  data_envi_2 = data_envi1[c("BO22_curvelmean_bdmean", "BO22_dissoxmean_bdmean",
                             "BO22_ironmean_bdmean" , "BO22_lightbotmean_bdmean",
                             "BO22_ph" , "BO22_nitratemean_bdmean" , "BO22_tempmean_bdmean",
                             "BO22_salinitymean_bdmean", "BO22_silicatemean_bdmean")]
  #Environmental distance
  dist.env = dist(data_envi_2, method = "euclidean")
  env = melt(as.matrix(dist.env))
  
  distances= data.frame(bc_dist,geo$value,env$value)
  colnames(distances) = c("Sample_1","Sample_2","BC_Dissimilarity","Distance_km","Distance_envi")
  
  #Log transformation
  distances$Distance_km_1= distances$Distance_km + 1
  distances$logDist=log10(distances$Distance_km_1)
  
  distances$Distance_envi_1 = distances$Distance_envi + 1
  distances$logEnv = log10(distances$Distance_envi_1)
  
  distances$Sample_1 = as.character(distances$Sample_1)
  distances$Sample_2 = as.character(distances$Sample_2)
  
  #Remove comparisons between same samples and "doubleton" comparisons
  distances2 = distances %>% dplyr::filter(Sample_1 < Sample_2)
  
  # Calculate slope
  lm_model_env <- lm(BC_Dissimilarity ~ logEnv, data = distances2)
  slopes_env[i] <- coef(lm_model_env)[2]
  
  lm_model_geo<- lm(BC_Dissimilarity ~ logDist, data = distances2)
  slopes_geo[i] <- coef(lm_model_geo)[2]
  
  #save in data frame
  slopes_data_frame = rbind(slopes_data_frame, data.frame(iterations=i,slope_ddr_env=slopes_env[i], slope_ddr_geo=slopes_geo[i]))
}

View(slopes_data_frame)
View(phylo_samples)

#Store the information in others vectors
#boostrap_slopes_sediment = slopes_data_frame
#boostrap_slopes_sediment$Habitat= c(rep("Sediment",1000))

#boostrap_slopes_content =slopes_data_frame
#boostrap_slopes_content$Habitat= c(rep("Gut content",1000))

#boostrap_slopes_tissue =slopes_data_frame
#boostrap_slopes_tissue$Habitat= c(rep("Gut tissue",1000))


#dataframe boostrap slopes
boostrap_slopes_df = rbind(boostrap_slopes_sediment,boostrap_slopes_content,boostrap_slopes_tissue)

mean(slopes_env)
mean(slopes_geo)

### DDR Geo
# mean sediment 0.1221501
# mean content 0.09975846
# mean tissue 0.04846488


### DDR Enviro
# mean sediment 0.497756
# mean content  0.4310629
# mean tissue 0.1911779

#DDR Geo slopes
ddr.slopes.geo = ggplot(boostrap_slopes_df,aes(x=Habitat, y=slope_ddr_geo, fill=Habitat)) +
  geom_point(aes(colour = Habitat), alpha=.7)+
  geom_boxplot(aes(fill = Habitat, colour=Habitat), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("") + 
  xlab("Habitat")+
  ylab("DDR slope")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut content"),
                                        c("Gut content","Gut tissue"),
                                        c("Sediment","Gut tissue")),
                      test="wilcox.test",
                      y_position = c(0.15,.16,0.17), 
                      map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.6) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment  ","  Gut content  ","  Gut tissue"))+
  scale_y_continuous(limits= c(0,0.2),breaks = c(0,0.1, 0.2))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black", vjust = 0.2),
        axis.title.x = element_text(size = 40, vjust = -1),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=38),
        legend.title = element_text(size=38),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
ddr.slopes.geo

#DDR Envi slopes
ddr.slopes.env = ggplot(boostrap_slopes_df,aes(x=Habitat, y=slope_ddr_env, fill=Habitat)) +
  geom_point(aes(colour = Habitat), alpha=.7)+
  geom_boxplot(aes(fill = Habitat, colour=Habitat), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("") + 
  xlab("Habitat")+
  ylab("DDR slope")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut content"),
                                        c("Gut content","Gut tissue"),
                                        c("Sediment","Gut tissue")),
                      test="wilcox.test",
                      y_position = c(0.6,.635,0.665), 
                      map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.6) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment  ","  Gut content  ","  Gut tissue"))+
  scale_y_continuous(limits= c(0,0.7),breaks = c(0,0.2,0.4, 0.6))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=48, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 38,color="black", vjust = 0.2),
        axis.title.x = element_text(size = 40, vjust = -1),
        axis.title.y = element_text(size=38),
        axis.text.y =element_text(size=38,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=38),
        legend.title = element_text(size=38),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
ddr.slopes.env


##### Panel Fig 4  ####
#remove legend form mantel tiss ddr to plot all together
plot_grid(mantel.plot.sed, mantel.plot.cont,mantel.plot.tiss,ddr.slopes.geo,
          mantel.plot.sed.env, mantel.plot.cont.env,mantel.plot.tiss.env,ddr.slopes.env, ncol=4, nrow=2)

#Save tiff 4000 2000

