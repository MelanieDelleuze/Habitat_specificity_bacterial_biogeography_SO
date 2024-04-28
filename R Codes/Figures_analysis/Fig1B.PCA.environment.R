##----------------------------------------##
##   Fig1B: PCA environmental variables   ##
##----------------------------------------##

# Set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")

# Loading required library
library(factoextra)
library(ggplot2)
library(ggtext)

# Environmental data fram: environmental conditions for each site 5 sites and 13 environmental variables from BIO-Oracle v.2.2
my_data_envi = read.table("envi.data.site.hab.spec.paper.txt",header=TRUE)
my_data_envi1 = my_data_envi[,-(1:6)] #remove columns that are not environmental data
# Rename rows with each Site names
rownames(my_data_envi1) = c("South Shetland I.", "Kerguelen I.", "Patagonia", "South Georgia", "Falkland/Malvinas I.")

# Perform a principal components analysis with all variables
pca.envi = prcomp(my_data_envi1, scale =TRUE)

# Testing correlation between environmental variables
core.data = cor(my_data_envi1, method = "pearson")
# Save correlation table 
write.table(core.data, file = "C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO/core.data.txt", sep = "\t",row.names = TRUE, col.names = NA)

# Leave only some variables removing positive correlated ones >0.95 to show only one vector for these variables in the PCA
my_data_envi2 = my_data_envi1[,-(10:11)]
my_data_envi2 = my_data_envi2[,-(5)]
my_data_envi2 = my_data_envi2[,-(11)]
my_data_envi2 = my_data_envi2[,-(1)]

# Rename variables cause some were almost same vector (see fviz_pca_biplot with pca.envi)
colnames(my_data_envi2) = c("Current velocity", "Oxygen", "Iron","Light" ,"pH",
                            "Nitrate/Phosphate", #Nitrate and Phosphote positive correlated >0.95 -> represented by the same arrow
                            "Temperature/Chlorophyl a*", #Temperature, Chlorophyl a, Primary Production, Carbon phytoplankton biomass,  positively correlated >0.95 -> represented by the same arrow
                            "Salinity","Silicate" )

# Perform a principal components analysis with these variables
pca.envi2 = prcomp(my_data_envi2, scale =TRUE)

# Plot PCA
pca.envi.plot = fviz_pca_biplot(pca.envi2, col.var = "cos2",label="var",labelsize = 22,
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                repel = TRUE, 
                                pointsize = 8, arrowsize = 4, alpha.var=0.7)+
  ggtitle("")+  labs(x = "PC1 (74%)", y = "PC2 (17%)")+
  geom_text_repel(aes(label = name), #limit overlap between labels
                  box.padding = 1, 
                  point.padding = 0.75, 
                  size = 22) +
  guides(colour = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme_minimal()+
  labs(colour=as.expression(bquote(bold("Contribution ("~cos^2~")")),size=22))+ #legend title
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x =  element_line(color="black"),
        axis.line.y =  element_line(color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 60,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=60,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y =element_blank(),
        legend.position="bottom",
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(2, 'cm'), 
        legend.key.width = unit(5, 'cm'),
        legend.text = element_text(size=58),
        legend.title = element_text(size=60, face="bold", hjust = 0.5),
        plot.margin = unit(c(2.6, 2.7, 0.5, 0.5),"cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
pca.envi.plot

# Export in tiff 3000 3100
