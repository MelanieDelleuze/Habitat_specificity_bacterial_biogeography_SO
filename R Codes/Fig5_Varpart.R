####-----------------------------##
#### Fig.5: Varpart analysis   ####
##-------------------------------##

#Load packages
library(phyloseq)
library(microViz)
library(geosphere)   
library(reshape) 
library(vegan)
library(ggplot2)
library(geosphere)
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

#### Geographic data ####
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

# Calculate the distances between sites in meters in a ellipsoid map WGS84
#(distGeo function does not provide a symmetric matrix as distm function) -> for loop to do it
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
#dist_geo_tiss=distance_matrix

#PCNM Principal Coordinates of Neighbourhood Matrix (Seppey et al., 2023)
#Sediment
geo.pcnm.sed = data.frame(scores(pcnm(dist_geo_sed)))
rownames(geo.pcnm.sed)= rownames(lats.sed2)
#Gut content
geo.pcnm.cont = data.frame(scores(pcnm(dist_geo_cont)))
rownames(geo.pcnm.cont)= rownames(lats.cont2)
#Gut tissue
geo.pcnm.tiss = data.frame(scores(pcnm(dist_geo_tiss)))
rownames(geo.pcnm.tiss)= rownames(lats.tiss2)

#### Environmental data ####
#Sediment
data_envi_sed = data.frame(sediment.rel@sam_data)
data_envi_sed = data_envi_sed[,-(1:15)]
#13 environmental variables
## PCA with all envi variables
pca.sed=prcomp(data_envi_sed, scale=TRUE)
## PCs per sample
PC.sed=data.frame(pca.sed$x)
# Creating a environmental distance matrix with PCA coordinates : PC1 & PC2
pc.scores.sed.pc = PC.sed[-(3:ncol(PC.sed))]

#Gut content
data_envi_cont = data.frame(content.rel@sam_data)
data_envi_cont = data_envi_cont[,-(1:15)]
#13 environmental variables
## PCA with all envi variables
pca.cont=prcomp(data_envi_cont, scale=TRUE)
## PCs per sample
PC.cont=data.frame(pca.cont$x)
# Creating a environmental distance matrix with PCA coordinates : PC1 & PC2
pc.scores.cont.pc = PC.cont[-(3:ncol(PC.cont))]

#Gut tissue
data_envi_tiss = data.frame(tissue.rel@sam_data)
data_envi_tiss = data_envi_tiss[,-(1:15)]
#13 environmental variables
## PCA with all envi variables
pca.tiss=prcomp(data_envi_tiss, scale=TRUE)
## PCs per sample
PC.tiss=data.frame(pca.tiss$x)
# Creating a environmental distance matrix with PCA coordinates : PC1 & PC2
pc.scores.tiss.pc = PC.tiss[-(3:ncol(PC.tiss))]

#### Varpart analysis ####
#Sediment
geo.pcnm.sed2 = geo.pcnm.sed[c("PCNM1", "PCNM2")]
pc.scores.sed.pc2 = pc.scores.sed.pc[c("PC1", "PC2")]

var.sed = varpart(dist.sed,geo.pcnm.sed2, pc.scores.sed.pc2)

# Test statistically the results
#transform otu table (can't work with distance)
sed.otu= data.frame(t(sediment@otu_table))
# Testing Geography alone
anova.cca(rda(sed.otu,geo.pcnm.sed2,pc.scores.sed.pc2)) #***
# Testing Environment alone
anova.cca(rda(sed.otu,pc.scores.sed.pc2,geo.pcnm.sed2)) #***
#All fractions statistical significant

#Gut content
geo.pcnm.cont2= geo.pcnm.cont[c("PCNM1", "PCNM2")]
pc.scores.cont.pc2=pc.scores.cont.pc[c("PC1","PC2")]

var.cont = varpart(dist.cont,geo.pcnm.cont2, pc.scores.cont.pc2)

# Test statistically the results
#transform otu table (can't work with distance)
cont.otu= data.frame(t(content.rel@otu_table))
# Testing Geography alone
anova.cca(rda(cont.otu,geo.pcnm.cont2, pc.scores.cont.pc2)) #***
# Testing Environment alone
anova.cca(rda(cont.otu,pc.scores.cont.pc2,geo.pcnm.cont2)) #***
#All fractions statistical significant

#Gut tissue
geo.pcnm.tiss2= geo.pcnm.tiss[c("PCNM1", "PCNM2")]
pc.scores.tiss.pc2=pc.scores.tiss.pc[c("PC1","PC2")]

var.tiss = varpart(dist.tiss,geo.pcnm.tiss2, pc.scores.tiss.pc2)

# Test statistically the results
#transform otu table (can't work with distance)
tiss.otu= data.frame(t(tissue.rel@otu_table))
# Testing Geography alone
anova.cca(rda(tiss.otu,geo.pcnm.tiss2, pc.scores.tiss.pc2)) #***
# Testing Environment alone
anova.cca(rda(tiss.otu,pc.scores.tiss.pc2,geo.pcnm.tiss2)) #***
#All fractions statistical significant

#### Plot the results ####

#Extract the data and make data.frame
varpart.sed= data.frame(var.sed[["part"]]$indfract)
rownames(varpart.sed) = c("Geography", "Environment", "Geography & Environment", "Residuals")
varpart.sed = varpart.sed[-4,]
varpart.sed$Habitat = rep("Sediment" ,3)

varpart.cont= data.frame(var.cont[["part"]]$indfract)
rownames(varpart.cont) = c("Geography", "Environment", "Geography & Environment", "Residuals")
varpart.cont = varpart.cont[-4,]
varpart.cont$Habitat = rep("Gut content" ,3)

varpart.tiss= data.frame(var.tiss[["part"]]$indfract)
rownames(varpart.tiss) = c("Geography", "Environment", "Geography & Environment", "Residuals")
varpart.tiss = varpart.tiss[-4,]
varpart.tiss$Habitat = rep("Gut tissue" ,3)
#Global dataframe
varpart.df=  rbind(varpart.sed,varpart.cont,varpart.tiss)
varpart.df$Factor= rownames(varpart.df)
varpart.df$Factor = c("Geography" ,               "Environment" ,             "Geography & Environment", 
                      "Geography"  ,             "Environment" ,           
                      "Geography & Environment", "Geography" ,             
                      "Environment",             "Geography & Environment" )

#Adjust to have the values with 2 digits
varpart.df$Adj.R.squared = round(varpart.df$Adj.R.squared , digits = 2)
varpart.df$Significance = c("0.001", "0.001","1", "0.001", "0.001","1", "0.001", "0.001","1" )

#To order the plot 
varpart.df$Habitat <- factor(varpart.df$Habitat , levels=c("Sediment", "Gut content", "Gut tissue"))
varpart.df$Factor <- factor(varpart.df$Factor , levels=c("Geography", "Environment", "Geography & Environment"))

varpartplot = ggplot(varpart.df, aes(fill=Factor, y=Adj.R.squared, x=Habitat, label=Adj.R.squared)) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(size = 11, position = position_stack(vjust = 0.5))+
  annotate("text", x=1, y=0.69, label="***", size=11, hjust=-1)+
  annotate("text", x=1, y=0.23, label="***", size=11, hjust=-1)+
  annotate("text", x=2, y=0.535, label="***", size=11, hjust=-1)+
  annotate("text", x=2, y=0.295, label="***", size=11, hjust=-1)+
  annotate("text", x=3, y=0.21, label="***", size=11, hjust=-1)+
  annotate("text", x=3, y=0.075, label="***", size=11, hjust=-1)+
  scale_fill_manual(values = c("#4BAA7C", "#BC4E91", "#58B3C6"))+
  labs(fill="")+
  xlab("")+ylab("Explained variation (%)")+ggtitle("")+
  scale_y_continuous(limits= c(-0.1,1),breaks = c(0,0.5,1))+
  geom_hline(yintercept = 0,
             color = "black",
             linewidth=0.8,
             linetype = "dashed") +
  theme(panel.background = element_blank(),
        plot.title = element_text(size=42, hjust=0.5),
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'),
        legend.position=c(0.75,0.9),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))

varpartplot
#Save tiff 1500 1500
#pdf 15 15