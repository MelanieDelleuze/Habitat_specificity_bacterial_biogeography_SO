##-------------------------------------------------##
####    Figure 8: Figure ecological Processes    ####
##-------------------------------------------------##

library(phyloseq)
library(microViz)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggpubr)

#Set working directory
setwd("~/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#To calculate the processes for each habitat we follow the script developed by Schwob G. #LINK GITHUB

# Import the data of the ecological processes for each habitat intra and inter site
ecological_processes_habitat = read.table("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/ecological_processes_habitat.txt",header = T)
ecological_processes_habitat_counts = read.table("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/ecological_processes_habitat_counts.txt",header = T)

#### Within sites ####
intra_site= subset(ecological_processes_habitat, Dataset == "intra_site")

facet_order = c("Gut tissue","Gut content","Sediment")
intra_site$Type <- factor(intra_site$Type, levels = facet_order)

barplot_site = ggplot(intra_site, aes(y = mean_percent, 
                                              x = Type, 
                                              fill=factor(result, levels = c( "Ecological drift","Dispersal limitation","Homogenizing dispersal","Variable selection","Homogeneous selection")))) +
  geom_bar(stat = "identity", position = "stack", color="black") +
  scale_fill_manual(values = c("Homogenizing dispersal" = "#a8dadc",
                               "Dispersal limitation" = "#219ebc", 
                               "Ecological drift" = "#1d3557", 
                               "Variable selection" = "#E78C47", 
                               "Homogeneous selection" = "#f9c74f"),
                    breaks = c("Homogenizing dispersal","Homogeneous selection","Dispersal limitation", 
                               "Variable selection","Ecological drift"),
                    labels = c("Homogenizing dispersal","Homogeneous selection","Dispersal limitation", 
                               "Variable selection","Ecological drift"),
                    name=NULL)+
  #ylim(0, 100) +
  labs(y = "Ecological processes (%)") + ggtitle("A. Intra-Site")+ 
  theme_bw()+
  theme(plot.title = element_text(size=40, face="bold", color="black"),
        axis.title.x = element_text(size=30),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(angle = -90),
        legend.text=element_text(size=20),
        legend.position = 'bottom',
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.10, 'cm'),
        axis.text.x=element_text(size=30, color ="black"),
        #strip.background = element_blank(colour="black", fill="white", size=1.5, linetype="solid"),
        strip.text = element_text(size=30),
        axis.text.y=element_text(size=35, color ="black",face="bold"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill="NA", size=2))+
  coord_flip() +
  guides(fill=guide_legend(nrow=5, ncol = 2, byrow=TRUE)) +
  scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Habitat") #+
barplot_site

#### Among sites ####
inter_site= subset(ecological_processes_habitat, Dataset == "inter_site")

inter_site$Type <- factor(inter_site$Type, levels = facet_order)

barplot_inter_site = ggplot(inter_site, aes(y = mean_percent, 
                                                    x = Type, 
                                                    fill=factor(result, levels = c( "Ecological drift","Dispersal limitation","Homogenizing dispersal","Variable selection","Homogeneous selection")))) +
  geom_bar(stat = "identity", position = "stack", color="black") +
  scale_fill_manual(values = c("Homogenizing dispersal" = "#a8dadc",
                               "Dispersal limitation" = "#219ebc", 
                               "Ecological drift" = "#1d3557", 
                               "Variable selection" = "#E78C47", 
                               "Homogeneous selection" = "#f9c74f"),
                    breaks = c("Homogenizing dispersal","Homogeneous selection","Dispersal limitation", 
                               "Variable selection","Ecological drift"),
                    labels = c("Homogenizing dispersal","Homogeneous selection","Dispersal limitation", 
                               "Variable selection","Ecological drift"),
                    name=NULL)+
  #ylim(0, 100) +
  labs(y = "Ecological processes (%)") + ggtitle("B. Inter-Site")+ 
  theme_bw()+
  theme(plot.title = element_text(size=40, face="bold", color="black"),
        axis.title.x = element_text(size=35),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(angle = -90),
        legend.text=element_text(size=20),
        legend.position = 'bottom',
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.10, 'cm'),
        axis.text.x=element_text(size=30, color ="black"),
        strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"),
        strip.text = element_text(size=35),
        axis.text.y=element_blank(),
        axis.ticks.y =element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill="NA", size=2))+
  #plot.margin = unit(c(0,0,0,-0.8), "cm"))+
  #facet_grid(Type ~ ., scales='free_y', space = "free", margins = FALSE) +
  coord_flip() +
  guides(fill=guide_legend(nrow=5, ncol = 2, byrow=TRUE)) +
  scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Habitat") +
  theme(ggh4x.axis.nestline = element_line(size = 1))
barplot_inter_site

leg <- get_legend(barplot_inter_site)
leg2=as_ggplot(leg)

fig8<- plot_grid(
  barplot_site + theme(legend.position="none"),
  barplot_inter_site + theme(legend.position="none"),
  align = 'hv',
  ncol = 2, nrow = 1, vjust=1, label_size = 30)
fig8

#save tiff 3000 600

#### Statistical tests: results in Table S9 ####
# Intra-site
intra_site_counts= subset(ecological_processes_habitat_counts, Dataset == "intra_site")

#Change the ecological process and the habitat
a <- as.vector(subset(intra_site_counts, result == "Dispersal limitation" & Type == "Sediment", select = percent)$percent)
b <- as.vector(subset(intra_site_counts, result == "Dispersal limitation" & Type == "Gut tissue", select = percent)$percent)

# Perform Wilcoxon rank sum test
wilcox.test(a, b, alternative = "two.sided")


###### Test habitat inter######
inter_site_counts= subset(ecological_processes_habitat_counts, Dataset == "inter_site")

#select unique combinations 
#levels(as.factor(inter2.2$Conca))
vector =c("FALK_KER" ,"FALK_KGI" ,"FALK_PAT", "FALK_SOG" ,
          "KER_KGI",  "KER_PAT",  "KER_SOG",  
          "KGI_PAT",  "KGI_SOG",  
          "PAT_SOG")
inter_site_counts_2 = inter_site_counts[inter_site_counts$Conca %in% vector,]

#Change the ecological process and the habitat
c <- as.vector(subset(inter_site_counts_2, result == "Variable selection" & Type == "Gut tissue", select = percent)$percent)
d <- as.vector(subset(inter_site_counts_2, result == "Variable selection" & Type == "Gut content", select = percent)$percent)

# Perform Wilcoxon rank sum test
wilcox.test(c, d, alternative = "two.sided")
