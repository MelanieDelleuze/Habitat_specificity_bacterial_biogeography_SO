##---------------------------------##
####   Preparation of the data   ####
##---------------------------------##

#### 1. Load Data and make Phyloseq object based on taxonomy #### 

#set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")

# The imported OTU table was filtered by Bokulich filter (0.005%, Bokulich et al., 2013) and was rarefied to 6000 sequences per sample
library(phyloseq)
# Import otu_table and tax_table as matrix
otu = as.matrix(read.table("otu.table.habitat.specificity.txt", header = TRUE,row.names = 1))
tax = as.matrix(read.table("taxa.table.habitat.specificity.txt", header = TRUE,row.names = 1))
# Transform to otu_table, tax_table objects and create a phyloseq
otu.table = otu_table(otu, taxa_are_rows = TRUE)
tax.table = tax_table(tax)
# Create phyloseq object
physeq = phyloseq(otu.table, tax.table)

# Import metadata
# the order of the rows in the metadata must match the order of the columns in the OTU table
metadata = read.table("metadata.habitat.specificity.txt", header = TRUE,row.names = 1)
sampledata = sample_data(data.frame(metadata, row.names=sample_names(physeq), stringsAsFactors=FALSE))
# Import the tree
tree = ape::read.tree("tree.habitat.specificity.tree")

# Merge data and tree to create final phyloseq object
physeq_habitat_specificity = merge_phyloseq(physeq, sampledata, tree)
physeq_habitat_specificity

# Save phyloseq 
save(physeq_habitat_specificity, file = "physeq_habitat_specificity.RData")

#### 2. Load Data PICRUSt2 and make Phyloseq object #### 
# Import PICRUSt2 outputs
metacyc_tab = "path_abun_unstrat_descrip_metacyc_hab_spec.tsv" #output table with abundances of metacyc pathways from PICRUSt2 
sam_tab = "metadata_hab_spec_func.tsv" #same metadata as for taxo but saved in tsv format

#Read outputs from PICRUSt2 as phyloseq object
library(microbiomeMarker)
ps_metacyc = import_picrust2(metacyc_tab, sam_tab, trait = "PATHWAY")
# Indicate row and column names as sample names
rownames(ps_metacyc@sam_data) = ps_metacyc@sam_data$Group
colnames(ps_metacyc@otu_table)= rownames(ps_metacyc@sam_data)
# Transform into relative abundances 
ps_metacyc_hab_spec = microbiome::transform(ps_metacyc, "compositional")
# Save phyloseq
save(ps_metacyc_hab_spec, file = "ps_metacyc_hab_spec.RData")


#### 3. Creation of subsets ####
library(microViz)
# By habitat
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
content = ps_filter(physeq_habitat_specificity,Type=="Gut content")
tissue = ps_filter(physeq_habitat_specificity,Type=="Gut tissue")

#Predicted functional data
sediment_func = ps_filter(ps_metacyc_hab_spec,Type=="Sediment")
content_func = ps_filter(ps_metacyc_hab_spec,Type=="Gut content")
tissue_func = ps_filter(ps_metacyc_hab_spec,Type=="Gut tissue")

# By site
kgi = ps_filter(physeq_habitat_specificity, Site =="KGI") #Antarctica site
ker = ps_filter(physeq_habitat_specificity, Site =="KER") #Kerguelen site
pat = ps_filter(physeq_habitat_specificity, Site =="PAT") #Patagonia site
sog = ps_filter(physeq_habitat_specificity, Site =="SOG") #South Georgia site
falk= ps_filter(physeq_habitat_specificity, Site =="FALK") #Falkland/Malvinas site




