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

#### 2. Creation of subsets ####
library(microViz)
# By habitat
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
content = ps_filter(physeq_habitat_specificity,Type=="Gut content")
tissue = ps_filter(physeq_habitat_specificity,Type=="Gut tissue")

# By site
kgi = ps_filter(physeq_habitat_specificity, Site =="KGI") #Antarctica site
ker = ps_filter(physeq_habitat_specificity, Site =="KER") #Kerguelen site
pat = ps_filter(physeq_habitat_specificity, Site =="PAT") #Patagonia site
sog = ps_filter(physeq_habitat_specificity, Site =="SOG") #South Georgia site
falk= ps_filter(physeq_habitat_specificity, Site =="FALK") #Falkland/Malvinas site




