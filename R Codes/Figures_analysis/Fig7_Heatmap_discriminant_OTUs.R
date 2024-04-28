##------------------------------------------------##
####   Figure 7: Discriminant OTUs per habitat  ####
##------------------------------------------------##

#Load packages
library(microbiomeMarker)
library(phyloseq)
library(RColorBrewer)
#Set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#Load data
load("physeq_habitat_specificity.RData") 

### LDA: Discriminant OTUs per Habitat ####
#Do it once and save the data
#set.seed(123)
#lefse_all= microbiomeMarker::run_lefse(physeq_habitat_specificty, group="Type",
#                                       bootstrap_n=100, lda_cutoff=2,
#                                       multigrp_strat = TRUE, taxa_rank="none")
#file edition to get the 10 most enriched functions in each species 
#lefse_all_dat = lefse_all@marker_table %>% group_by(enrich_group) %>% top_n(10, ef_lda)
#lefse_all_dat_mark =microbiomeMarker:: marker_table(lefse_all_dat)
#lefse_all_dat_mark_df = lefse_all_dat_mark #saving the lefse data for later

#lefse_all_dat_mark_df_2 = data.frame(lefse_all_dat_mark_df)
#write.table(lefse_all_dat_mark_df_2, file = "lefse_all_dat_mark_df.txt", sep = "\t",row.names = TRUE, col.names = NA)

#load table LDA
lefse_all_dat_mark_df = read.table("C:/Users/melad/Dropbox/Fondecyt Abatus 2021/6A. Estadia Melanie Delleuze/10.Paper/Script R/New_version_02_2024/lefse_all_dat_mark_df.txt", header=TRUE)
lefse_all_dat_mark_df = data.frame(lefse_all_dat_mark_df) #Contains the LDA score (Table S7)
#creating a vector
lefse.otu.hab = lefse_all_dat_mark_df[,2]
discriminant.otu.hab = as.vector(lefse.otu.hab)
otu.vector.hab = discriminant.otu.hab

#prune phyloseq with only discriminant OTUs from 3 habitats
hab_lda = prune_taxa(otu.vector.hab,physeq_habitat_specificty)

#In order to plot heatmap without gut content samples
#otu.vector.hab2 is the vector containing discriminant OTU from sediment and gut tissue
otu.vector.hab2 = c("Otu000007", "Otu000015", "Otu000029", "Otu000026", "Otu000083", "Otu000051" ,"Otu000079" ,"Otu000065", "Otu000025","Otu000073",
                    "Otu000012", "Otu000028", "Otu000033" ,"Otu000095", "Otu000047", "Otu000055", "Otu000203","Otu000074", "Otu000066", "Otu000041")
#prune just for sediment and Gut tissue samples
hab_lda2 = prune_taxa(otu.vector.hab2,physeq_habitat_specificty)
hab_lda2 = subset_samples(hab_lda2, Type != "Gut content") #remove gut content samples

#### Plot Heatmap ####

hab.lda.df.2 = data.frame(hab_lda2@otu_table)
hab.lda.taxa.df.2 = data.frame(hab_lda2@tax_table)
#following steps to improve the pheatmap labels
# add a column with the best taxonomical affiliation
hab.lda.taxa.df.2$Best_taxo = c("g;Woeseia OTU7", "g;Dethiosulfatarculus OTU74","c;Clostridia OTU55","f;Cyclobacteriaceae OTU79",
                                "p;Latescibacterota OTU47", "p;Latescibacterota OTU83","f;Thermoanaerobaculaceae subgroup 23 OTU65", "f;Thermoanaerobaculaceae subgroup 23 OTU29",
                                "k;Bacteria OTU66","f;Spirochaetaceae OTU33", "g;Spirochaeta_2 OTU41","g;Lutibacter OTU203",
                                "g;Lutibacter OTU12","k;Bacteria OTU95","c;Gammaproteobacteria OTU15" ,"c;Gammaproteobacteria OTU25",
                                "g;Delftia OTU28", "o;BD7-8 OTU73" ,"o;B2M28 OTU26", "g;Woeseia OTU51")
#access host vector hab_lda2@sam_data$Host 
#change last ones associated to the sediment as N.A. cause no host-associated
vector.host= c("Agassizii" , "Cavernosus" ,"Cavernosus", "Cavernosus", "Cavernosus", "Cavernosus", "Cordatus" ,  "Cavernosus",
               "Agassizii",  "Cavernosus" ,"Cordatus",   "Cavernosus" ,"Agassizii",  "Cordatus" ,  "Cavernosus", "Cavernosus",
               "Cavernosus", "Cavernosus", "Agassizii",  "Cavernosus", "Agassizii",  "Cavernosus", "Agassizii",  "Cordatus"  ,
               "Agassizii",  "Cavernosus", "Agassizii",  "Cavernosus", "Agassizii",  "Cavernosus" ,"Cavernosus", "Agassizii" ,
               "Agassizii" , "Agassizii",  "Agassizii",  "Agassizii",  "Cavernosus", "Agassizii",  "Cordatus" ,  "Agassizii" ,
               "Agassizii" , "Agassizii" , "Cavernosus", "Agassizii",  "Cordatus" ,  "Agassizii",  "Cavernosus", "Agassizii" ,
               "Agassizii" , "Cavernosus" ,"Cavernosus", "Agassizii",  "Agassizii",  "Cavernosus", "Agassizii",  "Agassizii" ,
               "Cavernosus", "Cavernosus", "Cavernosus", "NA" ,"NA" ,"NA", "NA" ,"NA", "NA" ,"NA", "NA" ,"NA", "NA", "NA", "NA", "NA", "NA" ,"NA", "NA", "NA", "NA", "NA", "NA" ,"NA",
               "NA", "NA" ,"NA" ,"NA" ,"NA" ,"NA" ,"NA")

#Rename appropiately
vector.host= gsub("Agassizii", "A. agassizii",vector.host)
vector.host=gsub("Cavernosus", "A. cavernosus",vector.host)
vector.host=gsub("Cordatus", "A. cordatus",vector.host)
#Change host vector
hab_lda2@sam_data$Host=vector.host

rownames(hab.lda.df.2) = hab.lda.taxa.df.2$Best_taxo 

#To add colors in row and col of the heatmap
my_sample_col2 <- data.frame(Site=hab_lda2@sam_data$Region,Host=hab_lda2@sam_data$Host,Habitat = hab_lda2@sam_data$Type)
row.names(my_sample_col2) <- colnames(hab.lda.df.2)

my_sample_row2 <- data.frame(Phylum = hab.lda.taxa.df.2$Phylum)
row.names(my_sample_row2) <- rownames(hab.lda.df.2)

#Vector of colors
my_colour2 = list(Habitat= c("Sediment"="#644536","Gut tissue"="#d4aa7d"),
                  Host=c("NA"="grey","A. agassizii" ="#84CAE7",  "A. cavernosus"= "#32785A", "A. cordatus"=  "#BD4F6C"),
                  Site = c("Antarctica"="#7AB8B2", "Falkland"="#E78419",
                           "Kerguelen"="#865AAC", "Patagonia"="#9D0208" , "South Georgia"="#90A955"),
                  Phylum = c("Acidobacteriota"="#ffcc00","Bacteria_unclassified"="grey","Bacteroidota"="#ff9900", "Desulfobacterota"="#ff6600",
                             "Firmicutes"="#669900", "Latescibacterota"="#99cc33",#"Modulibacteria"="#006699","Planctomycetota"="#3399cc",
                             "Proteobacteria" ="#3399cc", "Spirochaetota" ="#990066"))
breaklist = seq(-5, 5, by=0.5)

# Plot pheatmap
#run one first time with hab.lda.df.2 then change order of samples
ph.map.hab= pheatmap::pheatmap(hab.lda.df.2.order,filename = "C:/Users/melad/Desktop/figures.paper.test/Figure7/Figure7.heatmap.tiff",
                               scale="row", 
                               annotation_col = my_sample_col2,annotation_row = my_sample_row2,
                               annotation_colors = my_colour2,
                               cluster_cols=F,cluster_rows = T, 
                               show_colnames = F, border_color = TRUE,
                               color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaklist)), 
                               breaks = breaklist, 
                               treeheight_row=100,treeheight_col = 200,
                               legend = T,
                               fontsize = 35,cellwidth = 30, cellheight = 40,
                               width =60, height=30
)

#Changing order of samples
sample_order2 = c(#Sediment samples
                  "S1_Chi","S1_Chi_W" ,"S5_Chi","S4_Chi","S3_Chi", "S2_Chi",
                  "S4_SOG_22","S5_SOG_22" ,"S6_SOG_22", "S1_SOG_22", "S2_SOG_22", "S3_SOG_22"  ,
                  "S1_Paf" , "S2_Paf" ,  "S3_Paf" , "S4_Paf" , "S5_Paf"   ,
                  "S1_FALK_22", "S2_FALK_22",    "S3_FALK_22",  "S4_FALK_22", "S5_FALK_22",              
                  "S1_Par1"  ,  "S1_Par2"  , "S2_Par1"       , "S2_Par2" , "S3_Par1"  , "S3_Par2"  , 
                  #Gut tissue samples
                  "G10_Chi"  ,    "G14_Chi" ,  "G15_Chi"  ,   "G1_Chi"   ,   "G2_Chi"  ,  "G3_Chi"  ,    "G4_Chi"    ,      "G4_Chi_Tr" ,
                  "G5_Chi"        , "G6_Chi"    ,      "G6_Chi_Tr"    , "G7_Chi"      ,    "G7_Chi_Tr","G8_Chi"        , "G8_Chi_Tr" ,
                  "G9_Chi"   ,       "G9_Chi_Tr" ,  
                  "G2_PCR28_SOG_22" ,"G2_SOG_22", "G3_SOG_22" , "G4_PCR28_SOG_22",  "G4_SOG_22"  , "G5_PCR28_SOG_22", 
                  "G5_SOG_22"      , "G6_PCR28_SOG_22",  "G6_SOG_22"     ,
                  "G13_Paf"  , "G14_Paf"  ,"G15_Paf"   , "G2_Paf" ,"G5_Paf" ,  "G6_Paf" ,   
                  "G10_FALK_22",     "G12_FALK_22", "G13_FALK_22",   "G14_FALK_22", "G1_FALK_22", "G2_FALK_22","G3_FALK_22",    
                  "G5_FALK_22",  "G6_FALK_22", "G7_FALK_22","G8_FALK_22" ,   "G9_FALK_22" ,
                  "G10c_Par1" ,  "G12b_Par1" , "G13c_Par1" ,"G14c_Par1" ,  "G15c_Par1" ,      "G16c_Par1" , "G17c_Par1"  ,     "G18c_Par1",        
                  "G2c_Par1"  ,  "G3a_Par1"       , "G3c_Par1"  ,  "G6c_Par1"  ,   "G7c_Par1" ,  "G9a_Par1"      , 
                  "G9c_Par1" )
hab.lda.df.2.order <- hab.lda.df.2[ , sample_order2]
