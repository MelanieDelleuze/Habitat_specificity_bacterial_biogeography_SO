##----------------------------------------------------------##
####   Figure S6: Gut tissue Discriminant and Core OTUs   ####
##----------------------------------------------------------##

#Load packages
library(phyloseq)
library(microViz)
library(dplyr)
library(VennDiagram)
library(packcircles)
library(tibble)
library(gridExtra)

#Set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#Load data
load("physeq_habitat_specificity.RData")
#Make subsets
tissue = ps_filter(physeq_habitat_specificity,Type=="Gut tissue")
tissue.meta = data.frame(sample_data(tissue))

#Renaming
tissue.meta.new <- tissue.meta %>%
  mutate(Host = replace(Host, Host == "Agassizii", "A. agassizii"),
        Host = replace(Host, Host == "Cavernosus", "A. cavernosus"),
        Host = replace(Host, Host == "Cordatus", "A. cordatus") )

sample_data(tissue) = tissue.meta.new

#### Fig. S6A: Venn diagram gut tissue ####
png(file="C:/Users/melad/Desktop/figures.paper.test/FigureS6/Fig6A.venn.tiff",width=3000, height=3000, res=300)
MicEco::ps_venn(tissue,group = "Host", plot=TRUE, fraction=0,labels = list(cex = 2),
                   fills = list(fill=c("#E78419","#865AAC","#7AB8B2"),alpha=0.6),
                   quantities = list(type=c("percent","counts"), font = 1, cex=2),
                   adjust_labels=TRUE)
dev.off()

#### Fig.S6B :Core OTUs buble chart ####
G.tissue.core = core(tissue, detection = 0.1/100, prevalence = 65/100) #prev 85% 6 OTUs / prev 80% 10 OTUs/ prev 75% 12 OTUs / prev 70% 17 OTUs/prev 65% 28 OTUs
tissue.core = data.frame(G.tissue.core@tax_table)
write.table(tissue.core, file = "taxa.tissue.core.txt", sep = "\t",row.names = TRUE, col.names = NA)

##Gut tissue
core.tiss = as.data.frame(G.tissue.core@otu_table)
core.tiss$OTU_abundance= rowSums(core.tiss[1:nrow(core.tiss),])
core.tiss = cbind(rownames(core.tiss),core.tiss)
colnames(core.tiss)[1] <- "OTU"
core.tiss = core.tiss[,-(2:60)]
#rownames(core.tiss)
core.tiss$OTU.name=c( "OTU48", "OTU20", "OTU32", "OTU6", "OTU127", "OTU11", "OTU1", "OTU9", "OTU24",
                      "OTU67", "OTU85", "OTU4", "OTU14", "OTU18", "OTU12", "OTU5" ,"OTU19", "OTU22",
                      "OTU17", "OTU43", "OTU10", "OTU2", "OTU49", "OTU16", "OTU87", "OTU25" ,"OTU3",
                      "OTU26")

taxa.core.tiss= data.frame(G.tissue.core@tax_table)
taxa.core.tiss$Phylum

core.tiss = merge(core.tiss,taxa.core.tiss,by = 'row.names', all = TRUE)

# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing2 <- circleProgressiveLayout(core.tiss$OTU_abundance, sizetype='area')

# We can add these packing information to the initial data frame
core.tiss <- cbind(core.tiss, packing2)

# Check that radius is proportional to value. 
#We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
#plot(core.sed$radius, core.sed$OTU_abundance)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg2 <- circleLayoutVertices(packing2, npoints=50)
dat.gg2$phylum <- rep(core.tiss$Phylum, each=51)
dat.gg2$class <- rep(core.tiss$Class, each=51)

legend_title = "Phylum"

#Make the plot
circle.packing.chart.tiss = ggplot() + 
  # Make the bubbles
  geom_polygon(data = dat.gg2, aes(x, y, group = id, fill=phylum), colour = "black", alpha = 0.85) +
  scale_fill_manual(legend_title,
                    values = c("Actinobacteriota"="#1565c0",  "Bacteroidota"="#009688","Chloroflexi"="#2AB579",
                               "Desulfobacterota" ="#8bc34a",  "Myxococcota" ="#ff9800",
                               "Planctomycetota"="#f44336",   "Proteobacteria"="#ad1457",    "Spirochaetota"="#641E3C" ))+
  # Add text in the center of each bubble + control its size
  geom_text(data = core.tiss, aes(x, y, size=OTU_abundance, label = OTU.name)) +
  scale_size_continuous(range = c(3,9)) +
  # General theme:
  theme_void() + 
  guides(size="none")+theme(
    legend.text = element_text(size = 22, color="black"),
    legend.title = element_text(size = 22, color="black"),
    legend.key.size = unit(1, 'cm'), 
    legend.key.height = unit(0.8, 'cm'), 
    legend.key.width = unit(0.8, 'cm'))+
  coord_equal()


circle.packing.chart.tiss 
#Save tiff 

#### S6C : Heatmap LDA gut tissue ####
set.seed(123)
lefse_ab = microbiomeMarker::run_lefse(tissue, group="Host",
                                       bootstrap_n=100, lda_cutoff=2,
                                       multigrp_strat = TRUE, taxa_rank="none", kw_cutoff=0.01)
#file edition to get the 10 most enriched functions in each species 
lefse_ab_dat = lefse_ab@marker_table %>% group_by(enrich_group) %>% top_n(5, ef_lda)
lefse_ab_dat_mark = microbiomeMarker:: marker_table(lefse_ab_dat)
lefse_ab_dat_mark_df = lefse_ab_dat_mark #saving the lefse data for later
View(lefse_ab_dat_mark_df)

lefse_ab_dat_mark_df_2 = data.frame(lefse_ab_dat_mark_df) #saving the lefse data for later

write.table(lefse_ab_dat_mark_df_2, file = "lefse_ab_dat_mark_df_2.txt", sep = "\t",row.names = TRUE, col.names = NA)

#load data
lefse_ab_dat_mark_df = read.table("lefse_ab_dat_mark_df_2.txt", header=TRUE)
lefse_ab_dat_mark_df = data.frame(lefse_ab_dat_mark_df) #Contains the LDA score (Table S7)

#creating a vector
lefse.otu.ab = lefse_ab_dat_mark_df[,2]
discriminant.otu.ab = as.vector(lefse.otu.ab)
otu.vector.ab = discriminant.otu.ab

#prune sediment phyloseq with only discriminant OTUs
tiss_lda_host = prune_taxa(otu.vector.ab,tissue)

### Creating a table that contain the 10 most abundant OTU and the taxa information
lefse.otu.ab = lefse_ab_dat_mark[,1]
abundant.otu.ab = lefse.otu.ab$feature
otu.abund.ab = subset_taxa(tissue, rownames(tax_table(tissue)) %in% otu.vector.ab)
taxa.abund.ab = as.data.frame(otu.abund.ab@tax_table)
taxa.abund.ab = rownames_to_column(taxa.abund.ab,var="OTU")
taxa.abund.ab = taxa.abund.ab[match(otu.vector.ab, taxa.abund.ab$OTU),]
taxa.abund.ab$Site = c(rep("A. agassizii",5),rep("A. cavernosus",5),rep("A. cordatus",5))
View(taxa.abund.ab)

#Pheatmap function
tiss.lda.df = data.frame(tiss_lda_host@otu_table)
tiss.lda.taxa.df = data.frame(tiss_lda_host@tax_table)

tiss.lda.taxa.df$Best_taxo = c("f;Rhodobacteraceae OTU8", "c;Clostridia OTU690",
                               "p;Latescibacterota OTU47","k;Bacteria OTU66", "f;Spirochaetaceae OTU33",
                               "g;Spirochaeta_2 OTU4", "o;HOC36 OTU31", "g;Rhodopirellula OTU112", 
                               "g;Rhodopirellula OTU100","f;Pirellulaceae OTU61", 
                               "g;Pir4_lineage OTU16", "c;Gammaproteobacteria OTU44",
                               "g;Delftia OTU28", "g;Desulfoconvexum OTU30", "f;Desulfobacteraceae OTU3")


rownames(tiss.lda.df) = tiss.lda.taxa.df$Best_taxo 

my_sample_col <- data.frame(Host = tiss_lda_host@sam_data$Host)
row.names(my_sample_col) <- colnames(tiss.lda.df) 

my_sample_row <- data.frame(Phylum = tiss.lda.taxa.df$Phylum)
row.names(my_sample_row) <- rownames(tiss.lda.df)

my_colour = list(Host= c("A. agassizii"="#7AB8B2", "A. cavernosus"="#E78419",
                         "A. cordatus"="#865AAC"),
                 Phylum = c("Bacteria_unclassified"="grey",
                            "Desulfobacterota" ="#8bc34a",  "Firmicutes" ="#D5F05D","Latescibacterota"="#F7933C",
                            "Planctomycetota"="#f44336",   "Proteobacteria"="#ad1457",    "Spirochaetota"="#641E3C" ))

hm.gut =pheatmap::pheatmap(tiss.lda.df,filename = "C:/Users/melad/Desktop/figures.paper.test/FigureS6/Fig6C.tiff",
                           scale="row",annotation_col = my_sample_col,annotation_colors = my_colour,
                           annotation_row = my_sample_row,
                           show_colnames = F, border_color = TRUE,
                           color=colorRampPalette(c("blue", "white", "red"))(50),
                           treeheight_row=80,
                           treeheight_col =80,
                           legend = T,
                           fontsize = 22,cellwidth = 20, cellheight = 20,   width =30, height=10)

