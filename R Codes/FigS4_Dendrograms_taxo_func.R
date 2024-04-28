##------------------------------------------------------ ---##
#### Fig.S4 : Dendrogram Taxonomy-Functionality Sediment  ####
##----------------------------------------------------------##

# Load package
library(ggdendro)
library(dendextend)
library(phyloseq)
library(microViz)
# Load data
load("physeq_habitat_specificity.RData") 
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
sediment.rel = microbiome::transform(sediment,"compositional")
load("ps_metacyc_hab_spec.RData") 
sediment_func = ps_filter(ps_metacyc_hab_spec,Type=="Sediment")

#Sediment comparison bray-curtis taxonomy vs functionnality
dist.sed.mc = phyloseq::distance(sediment_func, method="bray")
dist.sed = phyloseq::distance(sediment.rel, method="bray")

dendro.sed <- as.dendrogram(hclust(dist.sed.mc), method="ward.D2")
dendro.sed
dendro.sed.plot <- ggdendrogram(dendro.sed)


# Make 2 dendrograms
d1 <- dist.sed.mc %>% dist() %>% hclust( method="ward.D") %>% as.dendrogram()
d2 <- dist.sed  %>% dist() %>% hclust( method="ward.D" ) %>% as.dendrogram()

# Custom these dendro, and place them in a list
dl <- dendlist(
  d2 %>% 
    set("labels_col", value = c("#7AB8B2","#9D0208","#90A955","#865AAC","#E78419"), k =5) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("#7AB8B2","#9D0208","#90A955","#865AAC","#E78419"), k =5),
  d1 %>% 
    set("labels_col", value = c("#90A955","#90A955","#9D0208","#7AB8B2","#E78419","#865AAC"), k =6) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("#90A955","#90A955","#9D0208","#7AB8B2","#E78419","#865AAC"), k =6))%>% 
  untangle(method = "step2side") # Find the best alignment layout


# Plot them together
tiff(file="C:/Users/melad/Desktop/figures.paper.test/FigureS4/Dendrograms_Sed_Taxo_Funct_bray.tiff",width=3500, height=2500, res=300)
tanglegram(dl, 
           common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=0.5,margin_outer=1,intersecting=TRUE,match_order_by_labels=TRUE, 
           color_lines = c( "#90A955","#90A955","#90A955","#90A955","#90A955","#90A955",
                            "#9D0208","#9D0208","#9D0208","#9D0208","#9D0208","#9D0208",
                            "#E78419","#E78419","#E78419","#E78419","#E78419",
                            "#865AAC","#865AAC","#865AAC","#865AAC","#865AAC",
                            "#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2","#7AB8B2"),
           lwd=2, edge.lwd=3.5,
           main_left= "Taxonomy", main_right = "Functionality", type="r", 
           cex.main=4,cex.axis=2.5,columns_width=c(10,2,10))

dev.off()
