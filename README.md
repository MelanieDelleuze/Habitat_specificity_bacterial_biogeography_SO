**Habitat_specificity_bacterial_biogeography_SO**
This is the associated R codes and functions used in the manuscript entitled "Habitat specificity modulates bacterial biogeographic patterns in the Southern Ocean", Delleuze et al. (2024) FEMS Microbiology Ecology.

The code was developped by Mélanie Delleuze and Guillaume Schwob.

The codes were divided for each figure and/or analyses: 

- **0.Extraction_environmental_data_Bio-oracle_v2.2.R**: Script to extract the environmental data from BIO-Oracle (https://www.bio-oracle.org/) based on the GPS coordinates of each sampling point.
- **1.Preparation_data.R** : Importation of tables and construction of phyloseq object and subsets.
- **Fig1B.PCA.environment.R** : Construction of the PCA based on environmental variables
- **Fig2_and_associated_supp_FigS1_S2_S3.R** : Alpha diversity, Phylogenetic Diversity, Levins' Niche Breadth index and betadisper analysis
- **Fig3_PCoA_Taxo_Func.R** : PCoA and PERMANOVA
- **Fig4_Varpart.R** : Variation partitionning analysis
- **Fig5_DDR.R** : Distance Decay Relationship
- **Fig6_Heatmap_discriminant_OTUs.R** : Linear discriminant analysis effect size (Lefse) by habitat and heatmap plotting
- **Fig7_Ecological_processes.R** : Script for the construction of the figure 7
- **Fig7_Ecological_processes_estimation.R** : Script for the estimation of the ecological processes developped by Guillaume Schwob and based on the scripts of Stegen et al. 2013 and Richter-Heitmann et al. 2020.
- **FigS4_Sediment_discriminant_OTU.R**: Lefse and core analysis of sediment communities
- **FigS5_Gut_tissue_discriminant_OTU.R** : Lefse and core analysis of gut tissue communities
- **FigS7_dendrograms_geo_envi.R** : Dendogramms based on geographical and environmental distances
