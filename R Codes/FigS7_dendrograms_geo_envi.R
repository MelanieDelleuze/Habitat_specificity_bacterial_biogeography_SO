##---------------------------------------------------------------##
####   Figure S7: Environmental and Geographical dendrograms   ####
##---------------------------------------------------------------##

#Load packages
library(ggdendro)
library(dendextend)
library(phyloseq)
library(microViz)
library(geosphere)

#Set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#Load data
load("physeq_habitat_specificity.RData")
#Make subsets
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
sediment.rel = microbiome::transform(sediment,"compositional")

#Bray-Curtis distance
dist.sed = phyloseq::distance(sediment.rel, method="bray")
#Geographical distance
long.sed = sediment.rel@sam_data$Longitude
lat.sed = sediment.rel@sam_data$Latitude
sample.sed = sediment.rel@sam_data$Group
lats.sed = data.frame(long.sed,lat.sed,sample.sed)
rownames(lats.sed) = lats.sed$sample.sed
lats.sed2= lats.sed[,-3]
#Distance with distGeo
data.coord =  lats.sed2   
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
dist_geo_sed=distance_matrix
#in km
dist_geo_sed_km = dist_geo_sed/1000

#Environmental distance
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

#environmental distance : dist.env.sed
#geographical distance : dist_geo_sed_km
#Bray-curtis distance : dist.sed

dist_geo_sed_km <- na.omit(dist_geo_sed_km)
dist_geo_sed_km =as.dist(dist_geo_sed_km)

# Make 2 dendrograms
d1.env <- dist.env.sed %>% dist() %>% hclust( method="ward.D") %>% as.dendrogram()
d2.geo <- dist_geo_sed_km  %>% dist() %>% hclust( method="ward.D" ) %>% as.dendrogram()
d3.bc <- dist.sed  %>% dist() %>% hclust( method="ward.D" ) %>% as.dendrogram()

# Custom these dendro, and place them in a list
dl2 <- dendlist(
  d2.geo %>% 
    set("labels_col", value = c("#865AAC","#90A955","#7AB8B2","#E78419","#9D0208"), k =5) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("#865AAC","#90A955","#7AB8B2","#E78419","#9D0208"), k =5),
  d3.bc %>% 
    set("labels_col", value = c("#7AB8B2","#9D0208","#90A955","#865AAC","#E78419"), k =5) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("#7AB8B2","#9D0208","#90A955","#865AAC","#E78419"), k =5))%>% 
  untangle(method = "step1") # Find the best alignment layout


#env "#9D0208","#7AB8B2","#90A955","#E78419","#865AAC"
#geo "#865AAC","#90A955","#7AB8B2","#E78419","#9D0208"
#bc  "#7AB8B2","#9D0208","#90A955","#865AAC","#E78419"

# Plot them together
tiff(file="C:/Users/melad/Desktop/figures.paper.test/FigureS7/Dendrograms_Sed_Geo_Taxo.tiff",width=3800, height=2500, res=300)
tanglegram(dl2, 
           common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=TRUE, 
           margin_inner=0.5,margin_outer=1,intersecting=TRUE,match_order_by_labels=TRUE, 
           color_lines = c(  "#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2",
                             "#9D0208","#9D0208","#9D0208","#9D0208","#9D0208","#9D0208",
                           "#E78419","#E78419","#E78419","#E78419","#E78419",
                          "#90A955","#90A955","#90A955","#90A955","#90A955","#90A955",
                           "#865AAC","#865AAC","#865AAC","#865AAC","#865AAC"),
           lwd=2, edge.lwd=3.5,
           main_left= "Geography", main_right = "Taxonomy", type="r", 
           cex.main=4,cex.axis=2.5,columns_width=c(12,4,12))

dev.off()

#Environment
#"#9D0208","#9D0208","#9D0208","#9D0208","#9D0208","#9D0208",
#"#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2",
#"#90A955","#90A955","#90A955","#90A955","#90A955","#90A955",
#"#865AAC","#865AAC","#865AAC","#865AAC","#865AAC",
#"#E78419","#E78419","#E78419","#E78419","#E78419"

#Geography
#"#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2",
#"#9D0208","#9D0208","#9D0208","#9D0208","#9D0208","#9D0208",
#"#E78419","#E78419","#E78419","#E78419","#E78419",
#"#90A955","#90A955","#90A955","#90A955","#90A955","#90A955",
#"#865AAC","#865AAC","#865AAC","#865AAC","#865AAC"