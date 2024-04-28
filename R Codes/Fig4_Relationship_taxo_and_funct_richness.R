##------------------------------------------------------------------------------------##
#### Fig.4: Relationship between the taxonomical and predicted functional richness  ####
##------------------------------------------------------------------------------------##

library(phyloseq)
library(metagMisc)
#Set working directory
setwd("C:/Users/melad/Documents/GitHub/Habitat_specificity_bacterial_biogeography_SO")
#Load data of taxonomy
load("physeq_habitat_specificity.RData")
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
content = ps_filter(physeq_habitat_specificity,Type=="Gut content")
tissue = ps_filter(physeq_habitat_specificity,Type=="Gut tissue")
#Load data from predicted fucntions
load("ps_metacyc_hab_spec.RData") 
#Creation subsets
sediment_func = ps_filter(ps_metacyc_hab_spec,Type=="Sediment")
content_func = ps_filter(ps_metacyc_hab_spec,Type=="Gut content")
tissue_func = ps_filter(ps_metacyc_hab_spec,Type=="Gut tissue")

#Sediment
#transform presence/absence to have number of OTUs per sample
sed_pa= phyloseq_standardize_otu_abundance(sediment, method="pa")
sed.otu.df= data.frame(t(sed_pa@otu_table))
sed.otu.df$Count_OTU = rowSums(sed.otu.df[,1:ncol(sed.otu.df)] )
sed.otu.df$Sample = rownames(sed.otu.df)
sed.otu.df =sed.otu.df[,-(1:1151)] #remove abundance data

#transform presence/absence to have number of OTUs per sample
sed_pa_mc= phyloseq_standardize_otu_abundance(sediment_func, method="pa")
sed.otu.df.mc= data.frame(t(sed_pa_mc@otu_table))
sed.otu.df.mc$Count_function = rowSums(sed.otu.df.mc[,1:ncol(sed.otu.df.mc)] )
sed.otu.df.mc$Sample_name = rownames(sed.otu.df.mc)
sed.otu.df.mc =sed.otu.df.mc[,-(1:382)] #remove abundance data

#Create data frame
sed.otu.mc.richness= merge(sed.otu.df,sed.otu.df.mc,by = 'row.names', all = TRUE)
rownames(sed.otu.mc.richness)=sed.otu.mc.richness$Row.names
sed.otu.mc.richness = sed.otu.mc.richness[,-1]
sed.otu.mc.richness = sed.otu.mc.richness[,-2]

#Fit linear model
fit.sed.rich <- lm(Count_OTU ~ Count_function  , data=sed.otu.mc.richness)
#view model summary
summary(fit.sed.rich)

# Plot
reg.rich.sed = ggplot(sed.otu.mc.richness, aes(y = Count_function, x=Count_OTU)) + 
  geom_point(size = 6, alpha = 0.8, colour = "#644536",shape = 21) + 
  geom_smooth(method = "lm", colour = "#644536", se=T,fill="#644536" ) +
  stat_regline_equation(label.y = 380, aes(label = ..eq.label..),size=11) +
  stat_cor(label.y = 370,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=11)+
  labs(y = "Functional richness", x = "Taxonomical richness") + 
  ggtitle("Sediment")+
  guides(fill="none")+
  scale_y_continuous(limits= c(280,380),breaks = c(280,330,380))+
  scale_x_continuous(limits= c(100,600),breaks = c(100,350,600))+
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

reg.rich.sed

#Repeat for other habitats
#Gut content 
#transform presence/absence to have number of OTUs per sample
cont_pa= phyloseq_standardize_otu_abundance(content, method="pa")
cont.otu.df= data.frame(t(cont_pa@otu_table))
cont.otu.df$Count_OTU = rowSums(cont.otu.df[,1:ncol(cont.otu.df)] )
cont.otu.df$Sample = rownames(cont.otu.df)
cont.otu.df =cont.otu.df[,-(1:1222)]

cont_pa_mc= phyloseq_standardize_otu_abundance(content_func, method="pa")
cont.otu.df.mc= data.frame(t(cont_pa_mc@otu_table))
cont.otu.df.mc$Count_function = rowSums(cont.otu.df.mc[,1:ncol(cont.otu.df.mc)] )
cont.otu.df.mc$Sample_name = rownames(cont.otu.df.mc)
cont.otu.df.mc =cont.otu.df.mc[,-(1:390)]

cont.otu.mc.richness= merge(cont.otu.df,cont.otu.df.mc,by = 'row.names', all = TRUE)
rownames(cont.otu.mc.richness)=cont.otu.mc.richness$Row.names
cont.otu.mc.richness = cont.otu.mc.richness[,-1]
cont.otu.mc.richness = cont.otu.mc.richness[,-2]

fit.cont.rich <- lm(Count_OTU ~ Count_function  , data=cont.otu.mc.richness)
#view model summary
summary(fit.cont.rich)


# Plot
reg.rich.cont = ggplot(cont.otu.mc.richness, aes(y = Count_function, x=Count_OTU)) + 
  geom_point(size = 6, alpha = 0.8, colour = "#9e6240",shape = 21) + 
  geom_smooth(method = "lm", colour = "#9e6240", se=T,fill="#9e6240" ) +
  stat_regline_equation(label.y = 380, aes(label = ..eq.label..),size=11) +
  stat_cor(label.y = 370,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=11)+
  labs(y = "Functional richness", x = "Taxonomical richness") + 
  ggtitle("Gut content")+
  guides(fill="none")+
  scale_y_continuous(limits= c(280,380),breaks = c(280,330,380))+
  scale_x_continuous(limits= c(100,600),breaks = c(100,350,600))+
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

reg.rich.cont

#Gut tissue
#transform presence/absence to have number of OTUs per sample
tiss_pa= phyloseq_standardize_otu_abundance(tissue, method="pa")
tiss.otu.df= data.frame(t(tiss_pa@otu_table))
tiss.otu.df$Count_OTU = rowSums(tiss.otu.df[,1:ncol(tiss.otu.df)] )
tiss.otu.df$Sample = rownames(tiss.otu.df)
tiss.otu.df =tiss.otu.df[,-(1:1172)]

tiss_pa_mc= phyloseq_standardize_otu_abundance(tissue_func, method="pa")
tiss.otu.df.mc= data.frame(t(tiss_pa_mc@otu_table))
tiss.otu.df.mc$Count_function = rowSums(tiss.otu.df.mc[,1:ncol(tiss.otu.df.mc)] )
tiss.otu.df.mc$Sample_name = rownames(tiss.otu.df.mc)
tiss.otu.df.mc =tiss.otu.df.mc[,-(1:392)]

tiss.otu.mc.richness= merge(tiss.otu.df,tiss.otu.df.mc,by = 'row.names', all = TRUE)
rownames(tiss.otu.mc.richness)=tiss.otu.mc.richness$Row.names
tiss.otu.mc.richness = tiss.otu.mc.richness[,-1]
tiss.otu.mc.richness = tiss.otu.mc.richness[,-2]


fit.tiss.rich <- lm(Count_OTU ~ Count_function  , data=tiss.otu.mc.richness)
#view model summary
summary(fit.tiss.rich)


# Plot
reg.rich.tiss = ggplot(tiss.otu.mc.richness, aes(y = Count_function, x=Count_OTU)) + 
  geom_point(size = 6, alpha = 0.8, colour = "#d4aa7d",shape = 21) + 
  geom_smooth(method = "lm", colour = "#d4aa7d", se=T,fill="#d4aa7d" ) +
  stat_regline_equation(label.y = 380, aes(label = ..eq.label..),size=11) +
  stat_cor(label.y = 370,aes(label = paste(..rr.label.. ,..p.label.., sep = "*`,`~")), size=11)+
  labs(y = "Functional richness", x = "Taxonomical richness") + 
  ggtitle("Gut tissue")+
  guides(fill="none")+
  scale_y_continuous(limits= c(280,380),breaks = c(280,330,380))+
  scale_x_continuous(limits= c(100,600),breaks = c(100,350,600))+
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

reg.rich.tiss

#Figure 4
Fig4 =plot_grid(reg.rich.sed , reg.rich.cont , reg.rich.tiss, ncol=3, nrow=1)
#Save tiff 3000 900
