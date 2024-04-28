##-----------------------------------------------------------##
#### Fig.3: PCoA along the gradient of habitat specificity ####
##-----------------------------------------------------------##

###### Fig.3A.PCoA per habitat Taxonomy ####
library(phyloseq)
library(microViz)
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

#Calculating distance Bray-Curtis
dist.sed = phyloseq::distance(sediment.rel, method="bray")
dist.cont = phyloseq::distance(content.rel, method="bray")
dist.tiss = phyloseq::distance(tissue.rel, method="bray")

#Ordinations
set.seed(123)
ordsed = ordinate(sediment.rel,method="PCoA", distance = dist.sed)
ordcont = ordinate(content.rel,method="PCoA", distance = dist.cont)
ordtiss = ordinate(tissue.rel,method="PCoA", distance = dist.tiss)

sed = plot_ordination(sediment.rel, ordsed, color= "Region",shape="Host") +
  ggtitle("Sediment")+
  scale_color_manual(breaks = c("Patagonia","Falkland","South Georgia","Kerguelen","Antarctica"),
                     values = c("#9D0208","#E78419","#90A955","#865AAC","#83C5BE"))+
  geom_point(size=6)+
  stat_ellipse(aes(group = Region), linetype = 2)+
  xlab("PCoA Axis 1 (31.5%)")+ #verify the axis contribution and modify accordingly
  ylab("PCoA Axis 2 (22.3%)")+
  scale_x_continuous(limits= c(-1,1),breaks = c(-1,0,1))+
  scale_y_continuous(limits= c(-0.6,0.6),breaks = c(-0.6,0,0.6))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=42, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))

sed

cont = plot_ordination(content.rel, ordcont, color= "Region", shape="Host") +
  ggtitle("Gut content")+
  scale_color_manual(breaks = c("Patagonia","Falkland","South Georgia","Kerguelen","Antarctica"),
                     values = c("#9D0208","#E78419","#90A955","#865AAC","#83C5BE"))+
  geom_point(size=6)+
  stat_ellipse(aes(group = Region), linetype = 2)+
  xlab("PCoA Axis 1 (34.0%)")+ #verify the axis contribution and modify accordingly
  ylab("PCoA Axis 2 (19.3%)")+
  scale_x_continuous(limits= c(-1,1),breaks = c(-1,0,1))+
  scale_y_continuous(limits= c(-0.6,0.6),breaks = c(-0.6,0,0.6))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=42, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
cont

tiss = plot_ordination(tissue.rel, ordtiss, color= "Region",shape="Host") +
  ggtitle("Gut tissue")+
  scale_color_manual(breaks = c("Patagonia","Falkland","Kerguelen","South Georgia","Antarctica"),
                     values = c("#9D0208","#E78419","#865AAC","#90A955","#83C5BE"),
                     labels=c("Patagonia","Falkland/Malvinas I.","Kerguelen I.","South Georgia","South Shetland I."))+
  scale_shape_discrete(labels=c(expression(italic("A.agassizii    ")),expression(italic("A.cavernosus")),expression(italic("A.cordatus     "))))+
  geom_point(size=6)+
  stat_ellipse(aes(group = Region), linetype = 2)+
  xlab("PCoA Axis 1 (20.7%)")+ #verify the axis contribution and modify accordingly
  ylab("PCoA Axis 2 (11.9%)")+
  labs(shape="Host",color="Site")+
  scale_x_continuous(limits= c(-1,1),breaks = c(-1,0,1))+
  scale_y_continuous(limits= c(-0.6,0.6),breaks = c(-0.6,0,0.6))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=42, hjust=0.5),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
tiss

#extracting the legend 
leg.pcoa <- get_legend(tiss)
leg.pcoa2=as_ggplot(leg.pcoa)

Pcoa = plot_grid(sed, cont,tiss, nrow=1, ncol=3)
#save svg tiff 3000 900

##### Fig.3B.PCoA per habitat Functions PICRUSt2 #####
#Load data
load("ps_metacyc_hab_spec.RData") 
#Creation subsets
sediment_func = ps_filter(ps_metacyc_hab_spec,Type=="Sediment")
content_func = ps_filter(ps_metacyc_hab_spec,Type=="Gut content")
tissue_func = ps_filter(ps_metacyc_hab_spec,Type=="Gut tissue")
# Transform in Relative abundances
sediment_func_rel = microbiome::transform(sediment_func,"compositional")
content_func_rel = microbiome::transform(content_func,"compositional")
tissue_func_rel = microbiome::transform(tissue_func,"compositional")

#Calculating distance BC
dist.sed.mc = phyloseq::distance(sediment_func_rel, method="bray")
dist.cont.mc = phyloseq::distance(content_func_rel, method="bray")
dist.tiss.mc = phyloseq::distance(tissue_func_rel, method="bray")

set.seed(123)
ordsed.mc = ordinate(sediment_func_rel,method="PCoA", distance = dist.sed.mc)
ordcont.mc = ordinate(content_func_rel,method="PCoA", distance = dist.cont.mc)
ordtiss.mc = ordinate(tissue_func_rel,method="PCoA", distance = dist.tiss.mc)

sed.mc = plot_ordination(sediment_func_rel, ordsed.mc, color= "Region",shape="Host") +
  ggtitle("")+
  scale_color_manual(breaks = c("Patagonia","Falkland","South Georgia","Kerguelen","Antarctica"),
                     values = c("#9D0208","#E78419","#90A955","#865AAC","#83C5BE"))+
  geom_point(size=6)+
  stat_ellipse(aes(group = Region), linetype = 2)+
  xlab("PCoA Axis 1 (42.6%)")+ #check and modify axis contribution accordingly
  ylab("PCoA Axis 2 (37.0%)")+
  scale_x_continuous(limits= c(-0.2,0.2),breaks = c(-0.2,0,0.2))+
  scale_y_continuous(limits= c(-0.2,0.2),breaks = c(-0.2,0,0.2))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=42),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size =34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
sed.mc

cont.mc = plot_ordination(content_func_rel, ordcont.mc, color= "Region", shape="Host") +
  ggtitle("")+
  scale_color_manual(breaks = c("Patagonia","Falkland","South Georgia","Kerguelen","Antarctica"),
                     values = c("#9D0208","#E78419","#90A955","#865AAC","#83C5BE"))+
  geom_point(size=6)+
  stat_ellipse(aes(group = Region), linetype = 2)+
  xlab("PCoA Axis 1 (43.7%)")+
  ylab("PCoA Axis 2 (22.8%)")+
  scale_x_continuous(limits= c(-0.2,0.2),breaks = c(-0.2,0,0.2))+
  scale_y_continuous(limits= c(-0.2,0.2),breaks = c(-0.2,0,0.2))+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=42),
        legend.position = "none",
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
cont.mc

tiss.mc = plot_ordination(tissue_func_rel, ordtiss.mc, color= "Region",shape="Host") +
  ggtitle("")+
  scale_color_manual(breaks = c("Patagonia","Falkland","Kerguelen","South Georgia","Antarctica"),
                     values = c("#9D0208","#E78419","#865AAC","#90A955","#83C5BE"),
                     labels=c("Patagonia","Falkland/Malvinas I.","Kerguelen I.","South Georgia","South Shetland I."))+
  scale_shape_discrete(labels=c(expression(italic("A.agassizii    ")),expression(italic("A.cavernosus")),expression(italic("A.cordatus     "))))+
  geom_point(size=6)+
  stat_ellipse(aes(group = Region), linetype = 2)+
  xlab("PCoA Axis 1 (40.0%)")+
  ylab("PCoA Axis 2 (16.9%)")+
  labs(shape="Host",color="Site")+
  scale_x_continuous(limits= c(-0.2,0.2),breaks = c(-0.2,0,0.2))+
  scale_y_continuous(limits= c(-0.2,0.2),breaks = c(-0.2,0,0.2))+
  theme(panel.background = element_blank(),
        legend.position = "none",
        plot.title=element_text(size=42),
        panel.grid.major=element_blank(),
        axis.line = element_line(linetype=1,color="grey20"),
        axis.text.x = element_text(size = 34,color="black"),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size=34),
        axis.text.y =element_text(size=34,color="black"),
        legend.key = element_blank(),
        legend.text = element_text(size=34),
        legend.title = element_text(size=34),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill="NA", size=3))
tiss.mc

Pcoa_func = plot_grid(sed.mc,cont.mc,tiss.mc, ncol=3, nrow=1)

Fig3 = plot_grid(sed,cont,tiss,
                 sed.mc,cont.mc,tiss.mc,
                 ncol=3, nrow=2)
#Save tiff 3000 1800

##---------------##
#### PERMANOVA ####
##---------------##
###### Permanova Taxonomy ######
set.seed(123)
otu= abundances(tissue.rel) #change phyloseq for each habitat
meta= meta(tissue.rel)
#Permanova test 
permanova= adonis2(t(otu) ~ Site,data=meta,permutations=999,method="bray", by = "terms")
#Saving results
write.table(permanova, file = "C:/Users/melad/Desktop/figures.paper.test/Figure3/permanova.bc.tiss.site.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#Permanova with unifrac dist matrix
set.seed(123)
meta.uf <- as(sample_data(tissue.rel), "data.frame")
dist.uf <- phyloseq::distance(tissue.rel, method = "unifrac")
dist.uf2 = t(dist.uf)
permanova.UF = adonis2(dist.uf2 ~ Province_1, data = meta.uf, perm=999,p.adjust="holm", by="terms")
#Saving results
write.table(permanova.UF, file = "C:/Users/melad/Desktop/figures.paper.test/Figure3/permanova.uf.tiss.prov.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

###### Pairwise Permanova Taxonomy Bray-Curtis ######
### Preparing data for pairwise permanova ###
otu.hab = t(otu_table(tissue.rel))
otu.hab = data.frame(otu.hab)
meta.hab = data.frame(sample_data(tissue.rel))
meta.hab = meta.hab[,-(9:26)]
hab.df = merge(meta.hab,otu.hab,by = 'row.names', all = TRUE)

### Pairwise permanova ###
set.seed(123)
pairwise.permanova = funfuns::pairwise.adonis(hab.df[,12:ncol(hab.df)],hab.df$Region,sim.method = "bray", p.adjust.m = "holm", permutations = 999)

#Saving results
write.table(pairwise.permanova, file = "C:/Users/melad/Desktop/figures.paper.test/Figure3/PERMANOVA/Taxonomy/pairwise.permanova.bc.tiss.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

###### Permanova Predicted Functions Picrust2 ######
set.seed(123)
#Sediment
otu.mc = abundances(tissue_func)
meta.mc = meta(tissue_func)
permanova.mc = adonis2(t(otu.mc)~Site,data=meta.mc,permutations=999,method="bray")
#Saving results
write.table(permanova.mc, file = "C:/Users/melad/Desktop/figures.paper.test/Figure3/PERMANOVA/Predicted functionality/permanova.bc.tiss.func.site.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


### Preparing data for pairwise permanova ###
otu.mc.2 = t(otu_table(tissue_func))
otu.mc.2 = data.frame(otu.mc.2)
meta.mc.2 = data.frame(sample_data(tissue_func))
meta.mc.2 = meta.mc.2[,-(9:27)]
otu.hab.df = merge(meta.mc.2,otu.mc.2,by = 'row.names', all = TRUE)

### Pairwise permanova ###
set.seed(123)
pairwise.permanova.mc = funfuns::pairwise.adonis(otu.hab.df[,10:ncol(otu.hab.df)],otu.hab.df$Region,sim.method = "bray", p.adjust.m = "holm", permutations = 999)
#Saving results
write.table(pairwise.permanova.mc, file = "C:/Users/melad/Desktop/figures.paper.test/Figure3/PERMANOVA/Predicted functionality/permanova.pairwise.bc.tiss.func.site.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


##----------------------------------------------------------##
#### Table S5: Mean Bray-Curtis differences between sites ####
##----------------------------------------------------------##

# Calculate distances
m = "bray"
p = content_func # change p by the phyloseq (sediment / tissue /content)
bray = phyloseq::distance(p, m)
bray.matrix = melt(as.matrix(bray))

bray.matrix = bray.matrix %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
s = "Group"
d = "Region"
sd = phyloseq::sample_data(p) 
sd = data.frame(sd) %>%
  dplyr::select(s,d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Region1")
bc.diss = left_join(bray.matrix, sd, by = "Var1")

colnames(sd) = c("Var2", "Region2")
bc.diss = left_join(bc.diss, sd, by = "Var2")

#write.table(bc.diss, file = "bc.diss.sed.txt", sep = "\t",row.names = TRUE, col.names = NA)

#Calculating mean bray-curtis dissimilarity between regions
#verification même valeur car comparaison 2 a 2 des echantillons apparaissent 2 fois
ant.ker.m = mean(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="Kerguelen"])
ant.pat.m = mean(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="Patagonia"])
ant.sog.m = mean(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="South Georgia"])
ant.falk.m = mean(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="Falkland"])
ker.pat.m = mean(bc.diss$value[bc.diss$Region1=="Kerguelen" & bc.diss$Region2=="Patagonia"])
ker.sog.m = mean(bc.diss$value[bc.diss$Region1=="Kerguelen" & bc.diss$Region2=="South Georgia"])
ker.falk.m = mean(bc.diss$value[bc.diss$Region1=="Kerguelen" & bc.diss$Region2=="Falkland"])
pat.sog.m = mean(bc.diss$value[bc.diss$Region1=="Patagonia" & bc.diss$Region2=="South Georgia"])
pat.falk.m = mean(bc.diss$value[bc.diss$Region1=="Patagonia" & bc.diss$Region2=="Falkland"])
sog.falk.m = mean(bc.diss$value[bc.diss$Region1=="South Georgia" & bc.diss$Region2=="Falkland"])

#Standard deviation
ant.ker.sd = sd(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="Kerguelen"])
ant.pat.sd = sd(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="Patagonia"])
ant.sog.sd = sd(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="South Georgia"])
ant.falk.sd = sd(bc.diss$value[bc.diss$Region1=="Antarctica" & bc.diss$Region2=="Falkland"])
ker.pat.sd = sd(bc.diss$value[bc.diss$Region1=="Kerguelen" & bc.diss$Region2=="Patagonia"])
ker.sog.sd = sd(bc.diss$value[bc.diss$Region1=="Kerguelen" & bc.diss$Region2=="South Georgia"])
ker.falk.sd = sd(bc.diss$value[bc.diss$Region1=="Kerguelen" & bc.diss$Region2=="Falkland"])
pat.sog.sd = sd(bc.diss$value[bc.diss$Region1=="Patagonia" & bc.diss$Region2=="South Georgia"])
pat.falk.sd = sd(bc.diss$value[bc.diss$Region1=="Patagonia" & bc.diss$Region2=="Falkland"])
sog.falk.sd = sd(bc.diss$value[bc.diss$Region1=="South Georgia" & bc.diss$Region2=="Falkland"])

#Stroe results
comparisons = c("Antarctica Kerguelen", "Antarctica Patagonia", "Antarctica South Georgia", "Antarctica Falkland",
                "Kerguelen Patagonia", "Kerguelen South Georgia","Kerguelen Falkland","Patagonia South Georgia",
                "Patagonia Falkland", "South Georgia Falkland")

mean_bc = c(ant.ker.m,ant.pat.m, ant.sog.m, ant.falk.m ,ker.pat.m ,
            ker.sog.m ,ker.falk.m ,  pat.sog.m ,pat.falk.m, sog.falk.m)

sd_bc = c(ant.ker.sd,ant.pat.sd, ant.sog.sd, ant.falk.sd ,ker.pat.sd ,
            ker.sog.sd ,ker.falk.sd ,  pat.sog.sd ,pat.falk.sd, sog.falk.sd)

#Create data frame
mean_bc = data.frame(comparisons, mean_bc, sd_bc)
#Saving results
write.table(mean_bc, file = "C:/Users/melad/Desktop/figures.paper.test/Figure3/PERMANOVA/mean_bc_cont_func.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
