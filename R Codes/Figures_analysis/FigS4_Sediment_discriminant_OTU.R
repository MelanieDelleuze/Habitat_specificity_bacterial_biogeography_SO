##--------------------------------------------------------##
####   Figure S5: Sediment Discriminant and Core OTUs   ####
##--------------------------------------------------------##

library(phyloseq)
library(microViz)
library(VennDiagram)
library(packcircles)
library(tibble)
library(gridExtra)

#Set working directory
setwd("~/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#Load data
load("physeq_habitat_specificity.RData")
#Make subsets
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
sediment.meta=data.frame(sample_data(sediment))

#Renaming
sediment.meta.new <- sediment.meta %>%
  mutate(Region = replace(Region, Region == "Antarctica", "South Shetland I."),
         Region = replace(Region, Region == "Falkland", "Falkland/Malvinas I."),
        Region = replace(Region, Region == "Kerguelen", "Kerguelen I."))
sample_data(sediment) = sediment.meta.new

#### Fig.S5A: Venn diagram ####
png(file="venn_sed_percents.png",width=3000, height=3000, res=300)
ps_venn(sediment,group = "Region", plot=TRUE, fraction=0.0001,labels = list(cex = 1.5),
        fills = list(fill=c("#E78419","#865AAC","#9D0208","#90A955","#7AB8B2"),alpha=0.6),
        quantities = list(type=c("percent"), font = 1, cex=1.5),
        adjust_labels=TRUE)
dev.off()

#### Fig.S5B: Core OTUs buble chart ####
sediment.core = core(sediment, detection = 0.1/100, prevalence = 90/100) #53 taxa 
# preparing data
core.sed = as.data.frame(sediment.core@otu_table)
core.sed = cbind(core.sed, c(rowSums(core.sed[1:nrow(core.sed),], na.rm = TRUE)))
colnames(core.sed)[29] <- "OTU_abundance"
core.sed = cbind(rownames(core.sed),core.sed)
colnames(core.sed)[1] <- "OTU"
core.sed = core.sed[,-(2:29)]
#To remove extra 000
core.sed$OTU.name=c("OTU157", "OTU7", "OTU149", "OTU152", "OTU108" ,"OTU20", "OTU32", "OTU6", "OTU8",
                    "OTU11", "OTU70","OTU79" ,"OTU1", "OTU35", "OTU36", "OTU136" ,"OTU101", "OTU102",
                    "OTU139", "OTU9", "OTU220" ,"OTU145", "OTU24", "OTU137", "OTU188", "OTU83", "OTU80",
                    "OTU65", "OTU29", "OTU71", "OTU266", "OTU14", "OTU12" ,"OTU5", "OTU22","OTU62",
                    "OTU17", "OTU10", "OTU2" ,"OTU16" ,"OTU332", "OTU40", "OTU58", "OTU15", "OTU25",
                    "OTU44", "OTU72", "OTU81", "OTU174", "OTU73", "OTU26" ,"OTU51" ,"OTU99")

taxa.core= data.frame(sediment.core@tax_table)
taxa.core$Phylum

core.sed = merge(core.sed,taxa.core,by = 'row.names', all = TRUE)

# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(core.sed$OTU_abundance, sizetype='area')

# We can add these packing information to the initial data frame
core.sed <- cbind(core.sed, packing)

# Check that radius is proportional to value. 
#We don't want a linear relationship, since it is the AREA that must be proportional to the value
#plot(core.sed$radius, core.sed$OTU_abundance)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
dat.gg$phylum <- rep(core.sed$Phylum, each=51)
dat.gg$class <- rep(core.sed$Class, each=51)

legend_title = "Phylum"

#Make the plot
circle.packing.chart.sed = ggplot() + 
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=phylum), colour = "black", alpha = 0.9) +
  scale_fill_manual(legend_title,
                    values = c("Acidobacteriota" = "#448aff",  "Actinobacteriota"="#1565c0",  "Bacteroidota"="#009688",
                               "Desulfobacterota" ="#8bc34a", "Latescibacterota" = "#ffc107", "Myxococcota" ="#ff9800",
                               "Planctomycetota"="#f44336",   "Proteobacteria"="#ad1457",    "Spirochaetota"="#641E3C" ),
                    labels=c("Acidobacteriota" , "Actinobacteriota", "Bacteroidota" ,    "Desulfobacterota", "Latescibacterota",
                             "Myxococcota",      "Planctomycetota",  "Proteobacteria",   "Spirochaetota" )) +
  # Add text in the center of each bubble + control its size
  geom_text(data = core.sed, aes(x, y, size=OTU_abundance, label = OTU.name)) +
  scale_size_continuous(range = c(1.5,8)) +
  # General theme:
  theme_void() + 
  guides(size="none")+theme(
    legend.text = element_text(size = 22, color="black"),
    legend.title = element_text(size = 22, color="black"),
    legend.key.size = unit(1, 'cm'), 
    legend.key.height = unit(0.8, 'cm'), 
    legend.key.width = unit(0.8, 'cm'))+
  coord_equal()


circle.packing.chart.sed 
#Save tiff 1000 1000

#### Fig.S5C: Heatmap LDA ####
#LDA Sediment
set.seed(123)
lefse_sed = microbiomeMarker::run_lefse(sediment, group="Region",
                                        bootstrap_n=100, lda_cutoff=2,
                                        multigrp_strat = TRUE, taxa_rank="none")
#file edition to get the 10 most enriched functions in each species 
lefse_sed_dat = lefse_sed@marker_table %>% group_by(enrich_group) %>% top_n(5, ef_lda)
lefse_sed_dat_mark =microbiomeMarker:: marker_table(lefse_sed_dat)
lefse_sed_dat_mark_df = lefse_sed_dat_mark #saving the lefse data for later
lefse_sed_dat_mark_df_2 = data.frame(lefse_sed_dat_mark) #saving the lefse data for later

write.table(lefse_sed_dat_mark_df_2, file = "lefse_sed_dat_mark_df.txt", sep = "\t",row.names = TRUE, col.names = NA)

#load data
lefse_sed_dat_mark_df = read.table("lefse_sed_dat_mark_df.txt", header=TRUE)
lefse_sed_dat_mark_df = data.frame(lefse_sed_dat_mark_df) #Contains the LDA score (Table S7)

#creating a vector
lefse.otu.sed = lefse_sed_dat_mark_df[,2]
discriminant.otu.sed = as.vector(lefse.otu.sed)
otu.vector.sed = discriminant.otu.sed

#prune sediment phyloseq with only discriminant OTUs
sed_lda = prune_taxa(otu.vector.sed,sediment)

### Creating a table that contain the 10 most abundant OTU and the taxa information
lefse.otu.sed = lefse_sed_dat_mark[,1]
abundant.otu.sed = lefse.otu.sed$feature
otu.abund.sed = subset_taxa(sediment, rownames(tax_table(sediment)) %in% abundant.otu.sed)
taxa.abund.sed = as.data.frame(otu.abund.sed@tax_table)
taxa.abund.sed = rownames_to_column(taxa.abund.sed,var="OTU")
taxa.abund.sed = taxa.abund.sed[match(abundant.otu.sed, taxa.abund.sed$OTU),]
taxa.abund.sed$Site = c(rep("Falkland/Malvinas I.",5),rep("Kerguelen I.",5),rep("Patagonia",5),rep("South Georgia",5), rep( "South Shetland I.",5))
View(taxa.abund.sed)
write.table(taxa.abund.sed, file = "taxa.abund.sed.txt", sep = "\t",row.names = TRUE, col.names = NA)

#Pheatmap function
sed.lda.df = data.frame(sed_lda@otu_table)
sed.lda.taxa.df = data.frame(sed_lda@tax_table)
#rownames(sed.lda.taxa.df) access rownmaes mand modify by best taxonomic affiliation
sed.lda.taxa.df$Best_taxo = c("g;Woeseia OTU7", "f;Sandaracinaceae OTU6","f;Sva1033 OTU121","g;Spirochaeta_2 OTU11",
                              "f;Cyclobacteriaceae OTU70", "f;Cyclobacteriaceae OTU79","g;Eudoraea OTU1",
                              "g;Maribacter OTU35", "g;Maribacter OTU36",
                              "g;Lutimonas OTU9","f;Bacteroidetes_BD2-2 OTU24",
                              "f;Pirellulaceae OTU67","f;Thermoanaerobaculaceae OTU29","g;Nitrospira OTU57","o:HOC36 OTU31",
                              "g;Blastopirellula OTU5","g;Blastopirellula OTU19","g;Bythopirellula OTU17","f;Pirellulaceae OTU10","g;Rubripirellula OTU2",
                              "c;Gammaproteobacteria OTU40","c;Gammaproteobacteria OTU15", "c;Gammaproteobacteria OTU25", 
                              "o;BD7-8 OTU73", "g;Woeseia OTU51")


rownames(sed.lda.df) = sed.lda.taxa.df$Best_taxo 

my_sample_col <- data.frame(Site = sed_lda@sam_data$Region)
row.names(my_sample_col) <- colnames(sed.lda.df)

my_sample_row <- data.frame(Phylum = sed.lda.taxa.df$Phylum)
row.names(my_sample_row) <- rownames(sed.lda.df)

my_colour = list(Site= c("South Shetland I."="#7AB8B2", "Falkland/Malvinas I."="#E78419",
                         "Kerguelen I."="#865AAC", "Patagonia"="#9D0208" , "South Georgia"="#90A955"),
                 Phylum = c("Acidobacteriota" = "#448aff",  "Bacteroidota"="#009688",
                            "Desulfobacterota" ="#8bc34a", "Myxococcota" ="#ff9800", "Nitrospirota"="#B03915",
                            "Planctomycetota"="#f44336",   "Proteobacteria"="#ad1457",    "Spirochaetota"="#641E3C" ))

#this pheatmap put antarctica between falk and ker 
ph.map.sed = pheatmap::pheatmap(sed.lda.df,#filename = "C:/Users/melad/Dropbox/Fondecyt Abatus 2021/6. Estadia Melanie Delleuze/10.Paper/Script R/Results/Figures.paper.09.2023.pdf/Heatmap.sediment.lda2.png",
                                scale="row",annotation_col = my_sample_col,annotation_colors = my_colour,
                                show_colnames = T, border_color = TRUE,
                                color=colorRampPalette(c("navy", "white", "red"))(50),
                                legend = T,
                                fontsize = 12,cellwidth = 20, cellheight = 20,
                                width =15, height=15)

#changing dendrogram order
col_dend <- ph.map.sed [[2]]
col_dend <- dendextend::rotate(col_dend, 
                               order = c("S3_SOG_22" , "S1_SOG_22","S2_SOG_22" ,"S6_SOG_22","S4_SOG_22","S5_SOG_22",
                                         "S3_FALK_22", "S4_FALK_22","S1_FALK_22","S2_FALK_22",
                                         "S5_FALK_22", 
                                         "S2_Par1" ,"S3_Par1","S1_Par2","S3_Par2","S1_Par1","S2_Par2" ,
                                         "S5_Paf","S1_Paf" ,"S3_Paf","S4_Paf" ,"S2_Paf", 
                                         "S1_Chi_W","S1_Chi","S2_Chi", "S3_Chi", "S4_Chi" ,"S5_Chi") )




# The pheatmap with the same clustering of heatmap
pheat.sed.clust = pheatmap::pheatmap(sed.lda.df, cluster_cols=as.hclust(col_dend),
                                     scale="row",annotation_col = my_sample_col,annotation_colors = my_colour,
                                     annotation_row = my_sample_row,
                                     show_colnames = F, border_color = TRUE,cutree_cols = 5,
                                     color=colorRampPalette(c("blue", "white", "red"))(50),
                                     treeheight_row=80,
                                     treeheight_col =80,
                                     legend = T,
                                     fontsize = 16,cellwidth = 20, cellheight = 20,width =20, height=10,
                                     filename = "~/GitHub/Habitat_specificity_bacterial_biogeography_SO/FigureS5/heatmap_lda_sed.tiff"
)
