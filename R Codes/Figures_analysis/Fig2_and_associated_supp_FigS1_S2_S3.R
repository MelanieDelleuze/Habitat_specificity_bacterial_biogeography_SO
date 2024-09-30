##----------------------------------------------------------------##
#### Fig.2: Characterization of the habitat specificity gradient ####
##----------------------------------------------------------------##

##### Fig.2A.Shannon alpha diversity #####

# Load required libraries
library(phyloseq)
#Set working directory
setwd("~/GitHub/Habitat_specificity_bacterial_biogeography_SO")
#Load phyloseq object
load("physeq_habitat_specificity.RData") 
physeq_habitat_specificity #in raw counts

# Estimate Shannon diversity index in each sample of the data set
richness.df = estimate_richness(physeq_habitat_specificity,measures=c("Shannon", "Simpson", "Chao1","Observed"))

# Create a table with the information of the samples and the measured richness index
richness.df2 = cbind(richness.df,sample_data(physeq_habitat_specificity))

# Test normality of the data
shapiro.test(richness.df2$Shannon)# p-values < 0.05, non normal distribution -> Non parametric statistics 
kruskal.test(richness.df2$Shannon, richness.df2$Type) #<2.2e-16

dunn_result = dunn.test::dunn.test(richness.df2$Shannon, richness.df2$Type,method="holm") #all <0.0000*

# Extract significant comparisons to add brackets to the boxplot
sig_comparisons <- as.data.frame(dunn_result)
sig_comparisons <- sig_comparisons[sig_comparisons$P.adjusted < 0.05, ]
sig_comparisons$Signifcance = c("****", "****", "****")

### Boxplot of alpha diversity ###
alpha.hab = ggplot(richness.df2,aes(x=Type, y=Shannon, fill=Type)) +
  geom_point(aes(colour = Type), alpha=.7)+
  geom_boxplot(aes(fill = Type, colour=Type), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  scale_y_continuous(breaks = c(1.5,3.5,5.5))+
  ggtitle("") + 
  xlab("")+
  ylab("Shannon diversity")+
  geom_signif(comparisons = list(c("Sediment", "Gut content"),
                                 c("Gut content", "Gut tissue"),
                                 c("Sediment", "Gut tissue")),
              annotations = sig_comparisons$Signifcance,
              y_position = c(5.3,5.5,5.7),
              map_signif_level = TRUE, vjust=0.7,textsize = 16)+
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    legend.key = element_blank(),#plot.margin = unit(c(1,1,1,1.5), "cm"),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
alpha.hab
#save tiff 1000 800


##### Fig.2B.Niche breadth index #####

### Niche breadth index for each habitat ###
library(MicroNiche)

# Creation of the subsets by habitat
sediment = ps_filter(physeq_habitat_specificity,Type=="Sediment")
content = ps_filter(physeq_habitat_specificity,Type=="Gut content")
tissue = ps_filter(physeq_habitat_specificity,Type=="Gut tissue")

# Calculating Bj index for each habitat separately
# Sediment
# Data frame in format needed by the function
data.otu.sed = data.frame(sediment@otu_table)
data.otu.sed1 = cbind(OTU=rownames(data.otu.sed),data.otu.sed)

# Vector with samples info
sampleInfo.sed = as.vector(sediment@sam_data$Region)

# Calculate Bj
Bj.sed = levins.Bn(data.otu.sed1,5,sampleInfo.sed)
# Remove OTUs below the limit of quantification (LOQ)
Bj.sed2 = subset(Bj.sed, Below.LOQ !="Y") #from 1151 OTUs to 201
# Mean value of Bj for sediment communities
Bcom.sed2 = mean(Bj.sed2$Bn) # 0.4191678

# Gut content
# Data frame in format needed by the function
data.otu.cont= data.frame(content@otu_table)
data.otu.cont1 = cbind(OTU=rownames(data.otu.cont),data.otu.cont)
# Vector with samples info
sampleInfo.cont =as.vector(content@sam_data$Region)
# Calculate Bj
Bj.cont = levins.Bn(data.otu.cont1,5,sampleInfo.cont)
# Remove OTUs below the limit of quantification (LOQ)
Bj.cont2 = subset(Bj.cont, Below.LOQ !="Y") #from 1222 OTUs to 200
# Mean value of Bj for gut content communities
Bcom.cont2 = mean(Bj.cont2$Bn) #0.3938988

# Gut tissue
# Data frame in format needed by the function
data.otu.tiss= data.frame(tissue@otu_table)
data.otu.tiss1 = cbind(OTU=rownames(data.otu.tiss),data.otu.tiss)
# Vector with samples info
sampleInfo.tiss = as.vector(tissue@sam_data$Region)
# Calculate Bj
Bj.tiss = levins.Bn(data.otu.tiss1,5,sampleInfo.tiss)
# Remove OTU below the limit of quantification (LOQ)
Bj.tiss2 = subset(Bj.tiss, Below.LOQ !="Y") #from 1172 OTUs to 154
# Mean value of Bj for gut content communities
Bcom.tiss2 = mean(Bj.tiss2$Bn) #0.3756253

### Making a dataframe with the results 
# Extracting Bn values
bn.sed = Bj.sed2$Bn 
bn.cont = Bj.cont2$Bn 
bn.tiss = Bj.tiss2$Bn 

#Extracting OTU names
otu.sed = rownames(Bj.sed2)
otu.cont = rownames(Bj.cont2)
otu.tiss = rownames(Bj.tiss2)

#Concatenate the vectors
otu.name = c(otu.sed,otu.cont,otu.tiss)
BN = c(bn.sed,bn.cont ,bn.tiss)
Habitat = c(rep("Sediment",201),rep("Gut content",200),rep("Gut tissue",154))

#Create a data frame with each Bn index per OTU in each compartment
levins.df = data.frame(cbind(otu.name,BN,Habitat))
levins.df$BN = as.numeric(levins.df$BN)
class(levins.df$BN)
#Save the dataframe
write.table(levins.df, file = "levins.niche.breadth.dataframe.new.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

### Statistic tests ###
# Shapiro-test
shapiro.test(levins.df$BN)
# p-value < 0.05, non normal distribution -> Non parametric statistics 

### Test stats Kruskall-Wallis + Dunn test Alpha diversity  ###
kruskal.test(levins.df,levins.df$BN~levins.df$Habitat) 
#Kruskal-Wallis rank sum test
#data:  levins.df
#Kruskal-Wallis chi-squared = 1142.3, df = 2, p-value < 2.2e-16

dunn_result_bn=dunn.test(levins.df$BN,levins.df$Habitat,method="holm")

# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_bn <- as.data.frame(dunn_result_bn)
sig_comparisons_bn <- sig_comparisons_bn[sig_comparisons_bn$P.adjusted < 0.05, ]
sig_comparisons_bn$Signifcance = c("**")


### Boxplot ###
levins.df$Habitat <- factor(levins.df$Habitat , levels=c("Sediment", "Gut content", "Gut tissue"))

# Boxplot Bj index
Levins.Bj = ggplot(levins.df,aes(x=Comp, y=BN, fill=Comp)) +
  geom_point(aes(colour = Comp), alpha=.7)+
  geom_boxplot(aes(fill = Comp, colour=Comp), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("") +
  xlab("")+
  ylab("Niche breadth index")+
  geom_signif(comparisons = list(c("Sediment","Gut tissue")),
              annotations = sig_comparisons_bn$Signifcance,
              y_position = c(0.96),
              map_signif_level = TRUE, vjust=0.7,textsize = 16)+
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  scale_y_continuous(breaks=c(0.2,0.6,1.0))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
Levins.Bj
#save tiff 1000 800

plot_grid(alpha.hab,Levins.Bj, ncol=1, nrow=2, labels="AUTO", label_size=40)
#save tiff 1000 1800


##-------------------------------------##
#### Fig.S1: Phylogenetic diversity: ####
##-------------------------------------##

##### Fig.S1.Faith's PD index #####
#BiocManager::install("genefilter")
#devtools::install_github("twbattaglia/btools")
library(btools)

meta.all = data.frame(physeq_habitat_specificity@sam_data)
#Estimate PD for each sample
all.pd = btools::estimate_pd(physeq_habitat_specificity)

# Data frame with PD index and information about the samples
pd_data = data.frame(all.pd, Type = meta.all$Type,Site=meta.all$Site, Region= meta.all$Region, Host = meta.all$Host)
View(pd_data)

# Statistical test
kruskal.test(pd_data$PD, pd_data$Type) #<2.2e-16

dunn_result_pd = dunn.test::dunn.test(pd_data$PD, pd_data$Type,method="holm") #all <0.0000*

# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_pd <- as.data.frame(dunn_result_pd)
sig_comparisons_pd <- sig_comparisons_pd[sig_comparisons_pd$P.adjusted < 0.05, ]
sig_comparisons_pd$Signifcance = c("****", "****", "****")


# Boxplot
pd.index = ggplot(pd_data,aes(x=Type, y=PD, fill=Type)) +
  geom_point(aes(colour = Type), alpha=.7)+
  geom_boxplot(aes(fill = Type, colour=Type), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("") +
  xlab("")+
  ylab("Phylogenetic Diversity")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut content"),
                                        c("Gut content","Gut tissue"),
                                        c("Sediment","Gut tissue")),
                      annotations = sig_comparisons_pd$Signifcance,
                      y_position = c(98,103,108), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_line(linetype=1,color="grey95"),
    axis.line = element_blank(),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
pd.index

#Save tiff 1000 800
#Save pdf 10 8 


##------------------------------------------------------------------------##
#### Fig.2: Characterization of the habitat specificity gradient by site ####
##------------------------------------------------------------------------##

##### Fig.S2A.Alpha diversity by site #####
# Subsets by sites
kgi = ps_filter(physeq_habitat_specificity, Site =="KGI")
ker = ps_filter(physeq_habitat_specificity, Site =="KER")
pat = ps_filter(physeq_habitat_specificity, Site =="PAT")
sog = ps_filter(physeq_habitat_specificity, Site =="SOG")
falk= ps_filter(physeq_habitat_specificity, Site =="FALK")

# Estimate Shannon diversity index in each sample of the data set and create a table with the data
rich.hab.site = estimate_richness(falk,measures="Shannon") #each time run it with different data set
sam.hab.site = cbind(rich.hab.site,sample_data(falk))

shapiro.test(sam.hab.site$Shannon)# p-values < 0.05, non normal distribution -> Non parametric statistics 
kruskal.test(sam.hab.site$Shannon, sam.hab.site$Type) #<2.2e-16

dunn_result_site = dunn.test(sam.hab.site$Shannon, sam.hab.site$Type,method="holm") #all different between them

# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_site <- as.data.frame(dunn_result_site)
sig_comparisons_site <- sig_comparisons_site[sig_comparisons_site$P.adjusted < 0.05, ]
sig_comparisons_site$Signifcance = c("**", "**", "****")

### Boxplots of alpha diversity by site ###

alpha.hab.falk= ggplot(sam.hab.site,aes(x=Type, y=Shannon, fill=Type)) +
  geom_point(aes(colour = Type), alpha=.7)+
  geom_boxplot(aes(fill = Type, colour=Type), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  scale_y_continuous(breaks = c(1.0,3.0,6.0))+
  ylim(1,6)+
  ggtitle("Falkland/Malvinas I.") + 
  xlab("")+
  ylab("")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut content"),
                                        c("Gut content","Gut tissue"),
                                        c("Sediment","Gut tissue")),
                      annotations = sig_comparisons_site$Signifcance,
                      y_position = c(5.3,5.5,5.7), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    plot.title = element_text(size=40),
    legend.key = element_blank(),#plot.margin = unit(c(1,1,1,1.5), "cm"),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
alpha.hab.falk
#save tiff 1000 800

FigS2A =plot_grid(alpha.hab.pat,alpha.hab.falk, alpha.hab.ker, 
          alpha.hab.sog,alpha.hab.kgi, 
          ncol=5, nrow=1, 
          #labels="AUTO", 
          label_size=40)
#save tiff 5000 1000


##### Fig.S2B. Niche breadth by site #####
### Calculating the Bj index in subset by site
###### South Shetland I. ######
kgi.sed = subset_samples(kgi,Type=="Sediment")
kgi.cont = subset_samples(kgi,Type=="Gut content")
kgi.tiss = subset_samples(kgi,Type=="Gut tissue")
#Sediment
data.otu.sed.kgi = data.frame(kgi.sed@otu_table)
data.otu.sed.kgi = cbind(OTU=rownames(data.otu.sed.kgi),data.otu.sed.kgi)
# Vector with sample info
sampleInfo.sed.kgi = as.vector(kgi.sed@sam_data$Group)
# Calculating Bj in sediment
Bj.sed.kgi = levins.Bn(data.otu.sed.kgi,6,sampleInfo.sed.kgi)
# Remove OTUs below the limit of quantification (LOQ)
Bj.sed.kgi2 = subset(Bj.sed.kgi, Below.LOQ !="Y") #from 1247 OTUs to 494
# Mean
Bcom.sed.kgi = mean(Bj.sed.kgi2$Bn) #0.5574483

# Content
data.otu.cont.kgi = data.frame(kgi.cont@otu_table)
data.otu.cont.kgi = cbind(OTU=rownames(data.otu.cont.kgi),data.otu.cont.kgi)
# vector with sample info
sampleInfo.cont.kgi = as.vector(kgi.cont@sam_data$Group)
Bj.cont.kgi = levins.Bn(data.otu.cont.kgi,36,sampleInfo.cont.kgi)
# Remove OTU below the limit of quantification (LOQ)
Bj.cont.kgi2 = subset(Bj.cont.kgi, Below.LOQ !="Y") #from 1247 OTUs to 652
# Mean
Bcom.cont.kgi = mean(Bj.cont.kgi2$Bn) #0.2681188

# Tissue
data.otu.tiss.kgi = data.frame(kgi.tiss@otu_table)
data.otu.tiss.kgi = cbind(OTU=rownames(data.otu.tiss.kgi),data.otu.tiss.kgi)
# vector with sample info
sampleInfo.tiss.kgi = as.vector(kgi.tiss@sam_data$Group)
Bj.tiss.kgi = levins.Bn(data.otu.tiss.kgi,17,sampleInfo.tiss.kgi)
# Remove OTU below the limit of quantification (LOQ)
Bj.tiss.kgi2 = subset(Bj.tiss.kgi, Below.LOQ !="Y") #from 1247 OTUs to 543
# Mean
Bcom.tiss.kgi = mean(Bj.tiss.kgi2$Bn) #0.1543377

### Making a dataframe with the resutls 
# Extracting Bn values
bn.sed.kgi = Bj.sed.kgi2$Bn
bn.cont.kgi = Bj.cont.kgi2$Bn
bn.tiss.kgi = Bj.tiss.kgi2$Bn

#Extracting OTU names
otu.sed.kgi = rownames(Bj.sed.kgi2)
otu.cont.kgi = rownames(Bj.cont.kgi2)
otu.tiss.kgi = rownames(Bj.tiss.kgi2)

#Concatenate the vectors
otu.name = c(otu.sed.kgi,otu.cont.kgi,otu.tiss.kgi)
BN = c(bn.sed.kgi, bn.cont.kgi, bn.tiss.kgi)
Comp = c(rep("Sediment",494),rep("Gut content",652),rep("Gut tissue",543))

#Create a data frame with each Bn index per OTU in each compartment
levins.df.kgi = data.frame(cbind(otu.name,BN,Comp))
levins.df.kgi$BN = as.numeric(levins.df.kgi$BN)
class(levins.df.kgi$BN)
#write.table(levins.df.kgi, file = "levinsBN.datafram.txt", sep = "\t",
#          row.names = TRUE, col.names = NA)


### Statistic tests ###
# Shapiro-test
shapiro.test(levins.df.kgi$BN)
# p-value < 0.05, non normal distribution -> Non parametric statistics 

### Test stats Kruskall-Wallis + Dunn test Alpha diversity  ###
kruskal.test(levins.df.kgi,levins.df.kgi$BN~levins.df.kgi$Comp) 
dunn_result_kgi = dunn.test(levins.df.kgi$BN,levins.df.kgi$Comp,method="holm")
# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_kgi <- as.data.frame(dunn_result_kgi)
sig_comparisons_kgi <- sig_comparisons_kgi[sig_comparisons_kgi$P.adjusted < 0.05, ]
sig_comparisons_kgi$Signifcance = c("****", "****", "****")

### Plot ###
levins.df.kgi$Comp <- factor(levins.df.kgi$Comp , levels=c("Sediment", "Gut content", "Gut tissue"))

# Boxplot Bj index
Levins.Bj.kgi = ggplot(levins.df.kgi,aes(x=Comp, y=BN, fill=Comp)) +
  geom_point(aes(colour = Comp), alpha=.7)+
  geom_boxplot(aes(fill = Comp, colour=Comp), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("South Shetland I.") +
  xlab("")+
  ylab("")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut tissue"),
                                        c("Sediment", "Gut content"),
                                        c("Gut content", "Gut tissue")),
                      annotations = sig_comparisons_kgi$Signifcance,
                      y_position = c(0.92, 0.96, 1), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  scale_y_continuous(breaks=c(0.5,1))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    plot.title = element_text(size=40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
Levins.Bj.kgi


###### Kerguelen I. ######
ker.sed = subset_samples(ker,Type=="Sediment")
ker.cont = subset_samples(ker,Type=="Gut content")
ker.tiss = subset_samples(ker,Type=="Gut tissue")
#Sediment
data.otu.sed.ker = data.frame(ker.sed@otu_table)
data.otu.sed.ker = cbind(OTU=rownames(data.otu.sed.ker),data.otu.sed.ker)
# vector with samples info
sampleInfo.sed.ker = as.vector(ker.sed@sam_data$Group)
# calculate niche breadth
Bj.sed.ker = levins.Bn(data.otu.sed.ker,5,sampleInfo.sed.ker)
#remove OTU below the limit of quantification (LOQ)
Bj.sed.ker2 = subset(Bj.sed.ker, Below.LOQ !="Y") #from 1247 OTUs to 616
# mean
Bcom.sed.ker = mean(Bj.sed.ker2$Bn) #0.513525

# Content
# Defining the filters
data.otu.cont.ker = data.frame(ker.cont@otu_table)
data.otu.cont.ker = cbind(OTU=rownames(data.otu.cont.ker),data.otu.cont.ker)
# vector
sampleInfo.cont.ker = as.vector(ker.cont@sam_data$Group)
#calculate bn
Bj.cont.ker = levins.Bn(data.otu.cont.ker,15,sampleInfo.cont.ker)
#remove OTU below the limit of quantification (LOQ)
Bj.cont.ker2 = subset(Bj.cont.ker, Below.LOQ !="Y") #from 1247 OTUs to 685
# Mean
Bcom.cont.ker = mean(Bj.cont.ker2$Bn) # 0.3633372

# Tissue
data.otu.tiss.ker = data.frame(ker.tiss@otu_table)
data.otu.tiss.ker = cbind(OTU=rownames(data.otu.tiss.ker),data.otu.tiss.ker)
# Vector with sample info
sampleInfo.tiss.ker = as.vector(ker.tiss@sam_data$Group)
# Calculating niche breadth
Bj.tiss.ker = levins.Bn(data.otu.tiss.ker,6,sampleInfo.tiss.ker)
#remove OTU below the limit of quantification (LOQ)
Bj.tiss.ker2 = subset(Bj.tiss.ker, Below.LOQ !="Y") #from 1247 OTUs to 378
# mean
Bcom.tiss.ker = mean(Bj.tiss.ker2$Bn) #0.2739333

### Making a dataframe with the resutls 
# Extracting Bn values
bn.sed.ker = Bj.sed.ker2$Bn
bn.cont.ker = Bj.cont.ker2$Bn
bn.tiss.ker = Bj.tiss.ker2$Bn

#Extracting OTU nams
otu.sed.ker = rownames(Bj.sed.ker2)
otu.cont.ker = rownames(Bj.cont.ker2)
otu.tiss.ker = rownames(Bj.tiss.ker2)

#Concatenate the vectors
otu.name = c(otu.sed.ker,otu.cont.ker,otu.tiss.ker)
BN = c(bn.sed.ker, bn.cont.ker, bn.tiss.ker)
Comp = c(rep("Sediment",616),rep("Gut content",685),rep("Gut tissue",378))

#Create a data frame with each Bn index per OTU in each compartment
levins.df.ker = data.frame(cbind(otu.name,BN,Comp))
levins.df.ker$BN = as.numeric(levins.df.ker$BN)
class(levins.df.ker$BN)
#write.table(levins.df.ker, file = "levinsBN.datafram.txt", sep = "\t",
#          row.names = TRUE, col.names = NA)

### Statistic tests ###
# Shapiro-test
shapiro.test(levins.df.ker$BN)
# p-value < 0.05, non normal distribution -> Non parametric statistics 

### Test stats Kruskall-Wallis + Dunn test Alpha diversity  ###
kruskal.test(levins.df.ker,levins.df.ker$BN~levins.df.ker$Comp) 
dunn_results_ker = dunn.test(levins.df.ker$BN,levins.df.ker$Comp,method="holm")
# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_ker <- as.data.frame(dunn_results_ker)
sig_comparisons_ker <- sig_comparisons_ker[sig_comparisons_ker$P.adjusted < 0.05, ]
sig_comparisons_ker$Signifcance = c("****", "****", "****")

### Plot ###
levins.df.ker$Comp <- factor(levins.df.ker$Comp , levels=c("Sediment", "Gut content", "Gut tissue"))

# Boxplot Bj index
Levins.Bj.ker = ggplot(levins.df.ker,aes(x=Comp, y=BN, fill=Comp)) +
  geom_point(aes(colour = Comp), alpha=.7)+
  geom_boxplot(aes(fill = Comp, colour=Comp), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("Kerguelen I.") +
  xlab("")+
  ylab("")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut tissue"),
                           c("Sediment", "Gut content"),
                           c("Gut content", "Gut tissue")),
                      annotations = sig_comparisons_ker$Signifcance,
                      y_position = c(0.92, 0.96, 1), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    plot.title = element_text(size = 40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
Levins.Bj.ker

###### Falkland/Malvinas I. ######
falk.sed = subset_samples(falk,Type=="Sediment")
falk.cont = subset_samples(falk,Type=="Gut content")
falk.tiss = subset_samples(falk,Type=="Gut tissue")
#Sediment
data.otu.sed.falk = data.frame(falk.sed@otu_table)
data.otu.sed.falk = cbind(OTU=rownames(data.otu.sed.falk),data.otu.sed.falk)
# vector
sampleInfo.sed.falk = as.vector(falk.sed@sam_data$Group)
# calculate niche breadth
Bj.sed.falk = levins.Bn(data.otu.sed.falk,5,sampleInfo.sed.falk)
#remove OTU below the limit of quantification (LOQ)
Bj.sed.falk2 = subset(Bj.sed.falk, Below.LOQ !="Y") #from 1247 OTUs to 741
# Mean
Bcom.sed.falk = mean(Bj.sed.falk2$Bn) #0.58239

# Content
data.otu.cont.falk = data.frame(falk.cont@otu_table)
data.otu.cont.falk = cbind(OTU=rownames(data.otu.cont.falk),data.otu.cont.falk)
# Vector of sample info
sampleInfo.cont.falk = as.vector(falk.cont@sam_data$Group)
# Calculating niche breadth
Bj.cont.falk = levins.Bn(data.otu.cont.falk,15,sampleInfo.cont.falk)
#remove OTU below the limit of quantification (LOQ)
Bj.cont.falk2 = subset(Bj.cont.falk, Below.LOQ !="Y") #from 1247 OTUs to 771
# Mean
Bcom.cont.falk = mean(Bj.cont.falk2$Bn) #0.3969634

# Tissue
data.otu.tiss.falk = data.frame(falk.tiss@otu_table)
data.otu.tiss.falk = cbind(OTU=rownames(data.otu.tiss.falk),data.otu.tiss.falk)
# Vector
sampleInfo.tiss.falk = as.vector(falk.tiss@sam_data$Group)
# calculating niche breadth
Bj.tiss.falk = levins.Bn(data.otu.tiss.falk,12,sampleInfo.tiss.falk)
#remove OTU below the limit of quantification (LOQ)
Bj.tiss.falk2 = subset(Bj.tiss.falk, Below.LOQ !="Y") #from 1247 OTUs to 648
# mean
Bcom.tiss.falk = mean(Bj.tiss.falk2$Bn) #0.2293684

### Making a dataframe with the resutls 
# Extracting Bn values
bn.sed.falk = Bj.sed.falk2$Bn
bn.cont.falk = Bj.cont.falk2$Bn
bn.tiss.falk = Bj.tiss.falk2$Bn

#Extracting OTU nams
otu.sed.falk = rownames(Bj.sed.falk2)
otu.cont.falk = rownames(Bj.cont.falk2)
otu.tiss.falk = rownames(Bj.tiss.falk2)

#Concatenate the vectors
otu.name = c(otu.sed.falk,otu.cont.falk,otu.tiss.falk)
BN = c(bn.sed.falk, bn.cont.falk, bn.tiss.falk)
Comp = c(rep("Sediment",741),rep("Gut content",771),rep("Gut tissue",648))

#Create a data frame with each Bn index per OTU in each compartment
levins.df.falk = data.frame(cbind(otu.name,BN,Comp))
levins.df.falk$BN = as.numeric(levins.df.falk$BN)
class(levins.df.falk$BN)

### Statistic tests ###
# Shapiro-test
shapiro.test(levins.df.falk$BN)
# p-value < 0.05, non normal distribution -> Non parametric statistics 

### Test stats Kruskall-Wallis + Dunn test Alpha diversity  ###
kruskal.test(levins.df.falk,levins.df.falk$BN~levins.df.falk$Comp) 
dunn_results_falk = dunn.test(levins.df.falk$BN,levins.df.falk$Comp,method="holm")
# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_falk <- as.data.frame(dunn_results_falk)
sig_comparisons_falk<- sig_comparisons_falk[sig_comparisons_falk$P.adjusted < 0.05, ]
sig_comparisons_falk$Signifcance = c("****", "****", "****")

### Plot ###
levins.df.falk$Comp <- factor(levins.df.falk$Comp , levels=c("Sediment", "Gut content", "Gut tissue"))

# Boxplot Bj index
Levins.Bj.falk = ggplot(levins.df.falk,aes(x=Comp, y=BN, fill=Comp)) +
  geom_point(aes(colour = Comp), alpha=.7)+
  geom_boxplot(aes(fill = Comp, colour=Comp), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("Falkland/Malvinas I.") +
  xlab("")+
  ylab("")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut tissue"),
                                        c("Sediment", "Gut content"),
                                        c("Gut content", "Gut tissue")),
                      annotations = sig_comparisons_kgi$Signifcance,
                      y_position = c(0.90, 0.94, 0.98), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    plot.title = element_text(size = 40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
Levins.Bj.falk

###### Patagonia #####
pat.sed = subset_samples(pat,Type=="Sediment")
pat.cont = subset_samples(pat,Type=="Gut content")
pat.tiss = subset_samples(pat,Type=="Gut tissue")
#Sediment
data.otu.sed.pat = data.frame(pat.sed@otu_table)
data.otu.sed.pat = cbind(OTU=rownames(data.otu.sed.pat),data.otu.sed.pat)
# Vector sample info
sampleInfo.sed.pat = as.vector(pat.sed@sam_data$Group)
# Calcualting niche breadth
Bj.sed.pat = levins.Bn(data.otu.sed.pat,6,sampleInfo.sed.pat)
#remove OTU below the limit of quantification (LOQ)
Bj.sed.pat2 = subset(Bj.sed.pat, Below.LOQ !="Y") #from 1247 OTUs to 819
# mean
Bcom.sed.pat = mean(Bj.sed.pat2$Bn) #0.5183589

# Content
data.otu.cont.pat = data.frame(pat.cont@otu_table)
data.otu.cont.pat = cbind(OTU=rownames(data.otu.cont.pat),data.otu.cont.pat)
# vector sample info
sampleInfo.cont.pat = as.vector(pat.cont@sam_data$Group)
# Calculating niche breadth
Bj.cont.pat = levins.Bn(data.otu.cont.pat,37,sampleInfo.cont.pat)
#remove OTU below the limit of quantification (LOQ)
Bj.cont.pat2 = subset(Bj.cont.pat, Below.LOQ !="Y") #from 1247 OTUs to 941
#mean
Bcom.cont.pat = mean(Bj.cont.pat2$Bn) #0.2366494

# Tissue
data.otu.tiss.pat = data.frame(pat.tiss@otu_table)
data.otu.tiss.pat = cbind(OTU=rownames(data.otu.tiss.pat),data.otu.tiss.pat)
# Vector sample info
sampleInfo.tiss.pat = as.vector(pat.tiss@sam_data$Group)
# Calculating niche breadth
Bj.tiss.pat = levins.Bn(data.otu.tiss.pat,15,sampleInfo.tiss.pat)
# Remove OTU below the limit of quantification (LOQ)
Bj.tiss.pat2 = subset(Bj.tiss.pat, Below.LOQ !="Y") #from 1247 OTUs to 718
# Mean
Bcom.tiss.pat = mean(Bj.tiss.pat2$Bn) # 0.1843909

### Making a dataframe with the resutls 
# Extracting Bn values
bn.sed.pat = Bj.sed.pat2$Bn
bn.cont.pat = Bj.cont.pat2$Bn
bn.tiss.pat = Bj.tiss.pat2$Bn

#Extracting OTU nams
otu.sed.pat = rownames(Bj.sed.pat2)
otu.cont.pat = rownames(Bj.cont.pat2)
otu.tiss.pat = rownames(Bj.tiss.pat2)

#Concatenate the vectors
otu.name = c(otu.sed.pat,otu.cont.pat,otu.tiss.pat)
BN = c(bn.sed.pat, bn.cont.pat, bn.tiss.pat)
Comp = c(rep("Sediment",819),rep("Gut content",941),rep("Gut tissue",718))

#Create a data frame with each Bn index per OTU in each compartment
levins.df.pat = data.frame(cbind(otu.name,BN,Comp))
levins.df.pat$BN = as.numeric(levins.df.pat$BN)
class(levins.df.pat$BN)

### Statistic tests ###
# Shapiro-test
shapiro.test(levins.df.pat$BN)
# p-value < 0.05, non normal distribution -> Non parametric statistics 

### Test stats Kruskall-Wallis + Dunn test Alpha diversity  ###
kruskal.test(levins.df.pat,levins.df.pat$BN~levins.df.pat$Comp) 
dunn_results_pat=dunn.test(levins.df.pat$BN,levins.df.pat$Comp,method="holm")
# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_pat <- as.data.frame(dunn_results_pat)
sig_comparisons_pat<- sig_comparisons_pat[sig_comparisons_pat$P.adjusted < 0.05, ]
sig_comparisons_pat$Signifcance = c("***", "****", "****") #arrange it manually 

### Plot ###
levins.df.pat$Comp <- factor(levins.df.pat$Comp , levels=c("Sediment", "Gut content", "Gut tissue"))
# Boxplot Bj index
Levins.Bj.pat = ggplot(levins.df.pat,aes(x=Comp, y=BN, fill=Comp)) +
  geom_point(aes(colour = Comp), alpha=.7)+
  geom_boxplot(aes(fill = Comp, colour=Comp), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("Patagonia") +
  xlab("")+
  ylab("Niche breadth index")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut tissue"),
                                        c("Sediment", "Gut content"),
                                        c("Gut content", "Gut tissue")),
                      annotations = sig_comparisons_pat$Signifcance,
                      y_position = c(0.90, 0.94, 0.98), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  scale_y_continuous(breaks=c(0.5,1))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    plot.title = element_text(size = 40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
Levins.Bj.pat

###### South Georgia #####
sog.sed = subset_samples(sog,Type=="Sediment")
sog.cont = subset_samples(sog,Type=="Gut content")
sog.tiss = subset_samples(sog,Type=="Gut tissue")
#Sediment
data.otu.sed.sog = data.frame(sog.sed@otu_table)
data.otu.sed.sog = cbind(OTU=rownames(data.otu.sed.sog),data.otu.sed.sog)
# Vector
sampleInfo.sed.sog = as.vector(sog.sed@sam_data$Group)
# Calculating niche breadth
Bj.sed.sog = levins.Bn(data.otu.sed.sog,6,sampleInfo.sed.sog)
#remove OTU below the limit of quantification (LOQ)
Bj.sed.sog2 = subset(Bj.sed.sog, Below.LOQ !="Y") #from 1247 OTUs to 785
# Mean
Bcom.sed.sog = mean(Bj.sed.sog2$Bn) #0.4412151

# Content
data.otu.cont.sog = data.frame(sog.cont@otu_table)
data.otu.cont.sog = cbind(OTU=rownames(data.otu.cont.sog),data.otu.cont.sog)
# vector
sampleInfo.cont.sog = as.vector(sog.cont@sam_data$Group)
# Calculating niche breadth
Bj.cont.sog = levins.Bn(data.otu.cont.sog,7,sampleInfo.cont.sog)
#remove OTU below the limit of quantification (LOQ)
Bj.cont.sog2 = subset(Bj.cont.sog, Below.LOQ !="Y") #from 1247 OTUs to 801
# Mean
Bcom.cont.sog = mean(Bj.cont.sog2$Bn) #0.3908323

# Tissue
data.otu.tiss.sog = data.frame(sog.tiss@otu_table)
data.otu.tiss.sog = cbind(OTU=rownames(data.otu.tiss.sog),data.otu.tiss.sog)
# vector
sampleInfo.tiss.sog = as.vector(sog.tiss@sam_data$Group)
# Calculating niche breadth
Bj.tiss.sog = levins.Bn(data.otu.tiss.sog,9,sampleInfo.tiss.sog)
#remove OTU below the limit of quantification (LOQ)
Bj.tiss.sog2 = subset(Bj.tiss.sog, Below.LOQ !="Y") #from 1247 OTUs to 560
# mean
Bcom.tiss.sog = mean(Bj.tiss.sog2$Bn)

### Making a dataframe with the resutls 
# Extracting Bn values
bn.sed.sog = Bj.sed.sog2$Bn
bn.cont.sog = Bj.cont.sog2$Bn
bn.tiss.sog = Bj.tiss.sog2$Bn

#Extracting OTU nams
otu.sed.sog = rownames(Bj.sed.sog2)
otu.cont.sog = rownames(Bj.cont.sog2)
otu.tiss.sog = rownames(Bj.tiss.sog2)

#Concatenate the vectors
otu.name = c(otu.sed.sog,otu.cont.sog,otu.tiss.sog)
BN = c(bn.sed.sog, bn.cont.sog, bn.tiss.sog)
Comp = c(rep("Sediment",785),rep("Gut content",801),rep("Gut tissue",560))

#Create a data frame with each Bn index per OTU in each compartment
levins.df.sog = data.frame(cbind(otu.name,BN,Comp))
levins.df.sog$BN = as.numeric(levins.df.sog$BN)
class(levins.df.sog$BN)


### Statistic tests ###
# Shapiro-test
shapiro.test(levins.df.sog$BN)
# p-value < 0.05, non normal distribution -> Non parametric statistics 

### Test stats Kruskall-Wallis + Dunn test Alpha diversity  ###
kruskal.test(levins.df.sog,levins.df.sog$BN~levins.df.sog$Comp) 
dunn_results_sog= dunn.test(levins.df.sog$BN,levins.df.sog$Comp,method="holm")
# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_sog<- as.data.frame(dunn_results_sog)
sig_comparisons_sog<- sig_comparisons_sog[sig_comparisons_sog$P.adjusted < 0.05, ]
sig_comparisons_sog$Signifcance = c("****", "****", "****")

### Plot ###
levins.df.sog$Comp <- factor(levins.df.sog$Comp , levels=c("Sediment", "Gut content", "Gut tissue"))

# Boxplot Bj index
Levins.Bj.sog = ggplot(levins.df.sog,aes(x=Comp, y=BN, fill=Comp)) +
  geom_point(aes(colour = Comp), alpha=.7)+
  geom_boxplot(aes(fill = Comp, colour=Comp), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("South Georgia") +
  xlab("")+
  ylab("")+
  ggpubr::geom_signif(comparison = list(c("Sediment","Gut tissue"),
                                        c("Sediment", "Gut content"),
                                        c("Gut content", "Gut tissue")),
                      annotations=sig_comparisons_sog$Signifcance,
                      y_position = c(0.92, 0.96, 1), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  ylim(0,1)+scale_y_continuous(breaks=c(0.5,1))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    plot.title = element_text(size = 40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
Levins.Bj.sog

FigS2B = plot_grid(Levins.Bj.pat,Levins.Bj.falk, Levins.Bj.ker, 
                   Levins.Bj.sog,Levins.Bj.kgi, 
                   ncol=5, nrow=1,
                   #labels="AUTO", 
                   label_size=40)


FigS2 = plot_grid(alpha.hab.pat,alpha.hab.falk, alpha.hab.ker, 
          alpha.hab.sog,alpha.hab.kgi,
          Levins.Bj.pat,Levins.Bj.falk, Levins.Bj.ker, 
          Levins.Bj.sog,Levins.Bj.kgi, 
          ncol=5, nrow=2, 
          #labels="AUTO", 
          label_size=40)
#save tiff 5000 1000


##---------------------------------------##
#### Fig.S3: Dispersion of the sammples ####
##---------------------------------------##

# Test point dispersion in space statistically
### Bray-curtis dissimilarity ###
dist.bc = phyloseq::distance(physeq_habitat_specificity, method="bray")
### Jaccard distance ###
dist.jc = phyloseq::distance(physeq_habitat_specificity, method="jaccard")
### Unifrac dissimilarity ###
dist.UF = phyloseq::distance(physeq_habitat_specificity, method="unifrac")

type = (physeq_habitat_specificity@sam_data)$Type
type = as.vector(type)
betadisp = betadisper(dist.bc,type)

shapiro.test(betadisp$distances)# p-values < 0.05, non normal distribution -> Non parametric statistics 
kruskal.test(betadisp$distances, betadisp$group) #<2.2e-16
dunn_resutls_betadisp = dunn.test(betadisp$distances, betadisp$group,method="holm") #all different between them
# Extract significant comparisons to add brackets to the boxplot
sig_comparisons_bd <- as.data.frame(dunn_resutls_betadisp)
sig_comparisons_bd<- sig_comparisons_bd[sig_comparisons_bd$P.adjusted < 0.05, ]
sig_comparisons_bd$Signifcance = c("****", "****")

#data.frame
df.betadisper = data.frame(dist = betadisp$distances,Type=betadisp$group)
df.betadisper$Type <- factor(df.betadisper$Type , levels=c("Sediment", "Gut content", "Gut tissue"))

#boxplot
betadisp.plot = ggplot(df.betadisper,aes(x=Type, y=dist, fill=Type)) +
  geom_point(aes(colour = Type), alpha=.7)+
  geom_boxplot(aes(fill = Type, colour=Type), alpha=.7) + 
  scale_color_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                     values = c("Sediment"="#644536","Gut content"= "#9e6240",
                                "Gut tissue"="#d4aa7d"))+
  scale_fill_manual(breaks = c("Sediment","Gut content","Gut tissue"),
                    values = c("Sediment"="#644536","Gut content"= "#9e6240",
                               "Gut tissue"="#d4aa7d"))+
  ggtitle("") +
  xlab("")+
  ylab("Variability")+
  ggpubr::geom_signif(comparison = list(c("Gut content","Gut tissue"),
                                        c("Sediment","Gut tissue")),
                      annotations=sig_comparisons_bd$Signifcance,
                      y_position = c(0.75, 0.79), map_signif_level=TRUE,
                      textsize=16, fontface = "italic", vjust=0.7) +
  scale_x_discrete(limits = c("Sediment","Gut content","Gut tissue"),
                   labels = c("Sediment","Gut content","Gut tissue"))+
  #scale_y_continuous(breaks = c(0,.5,1))+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid.major=element_blank(),
    axis.line = element_line(linetype=1,color="grey20"),
    axis.text.x = element_text(size=34,color="black"),
    axis.text.y = element_text(size=34,color="black"),
    axis.title.y = element_text(size=40),
    legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill="NA", size=3))
betadisp.plot

#save tiff 1000 900 and pdf 10 8

### For Table S3 ###
betadisp.bc = betadisper(dist.bc,type) #change type of distance
betadisp.jc = betadisper(dist.jc,type) #change type of distance
betadisp.uf = betadisper(dist.UF,type) #change type of distance

kruskal.test(betadisp.bc$distances, betadisp.bc$group) #<2.2e-16
dunn_resutls_betadisp_bc = dunn.test(betadisp.bc$distances, betadisp.bc$group,method="holm") #all different between them

kruskal.test(betadisp.jc$distances, betadisp.jc$group) #<2.2e-16
dunn_resutls_betadisp_jc = dunn.test(betadisp.jc$distances, betadisp.jc$group,method="holm") #all different between them

kruskal.test(betadisp.uf$distances, betadisp.uf$group) #<2.2e-16
dunn_resutls_betadisp_uf = dunn.test(betadisp.uf$distances, betadisp.uf$group,method="holm") #all different between them


#Mean
mean.sed.bc = mean(betadisp.bc$distances[betadisp.bc$group=="Sediment"])
mean.cont.bc = mean(betadisp.bc$distances[betadisp.bc$group=="Gut content"])
mean.tiss.bc = mean(betadisp.bc$distances[betadisp.bc$group=="Gut tissue"])
#Jaccard
mean.sed.jc = mean(betadisp.jc$distances[betadisp.jc$group=="Sediment"])
mean.cont.jc = mean(betadisp.jc$distances[betadisp.jc$group=="Gut content"])
mean.tiss.jc = mean(betadisp.jc$distances[betadisp.jc$group=="Gut tissue"])
#Unifrac
mean.sed.uf = mean(betadisp.uf$distances[betadisp.uf$group=="Sediment"])
mean.cont.uf = mean(betadisp.uf$distances[betadisp.uf$group=="Gut content"])
mean.tiss.uf = mean(betadisp.uf$distances[betadisp.uf$group=="Gut tissue"])

Habitat= c("Sediment", "Gut content", "Gut tissue","Sediment", "Gut content", "Gut tissue","Sediment", "Gut content", "Gut tissue")
Means = c(mean.sed.bc, mean.cont.bc, mean.tiss.bc,mean.sed.jc, mean.cont.jc, mean.tiss.jc,mean.sed.uf, mean.cont.uf, mean.tiss.uf)
comparisons = c(dunn_resutls_betadisp_bc$comparisons,dunn_resutls_betadisp_jc$comparisons,dunn_resutls_betadisp_uf$comparisons)
p_adjusted = c(dunn_resutls_betadisp_bc$P.adjusted,dunn_resutls_betadisp_jc$P.adjusted,dunn_resutls_betadisp_uf$P.adjusted)
Matrix = c(rep("Bray-Curtis",3),rep("Jaccard",3),rep("Unweighted UniFrac",3))


TableS3 = data.frame(Matrix, Habitat, Means, comparisons, p_adjusted)
write.table(TableS3, file = "~/GitHub/Habitat_specificity_bacterial_biogeography_SO/FigureS3/TableS3.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
