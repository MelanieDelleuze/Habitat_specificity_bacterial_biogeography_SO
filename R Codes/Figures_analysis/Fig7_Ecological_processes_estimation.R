##---------------------------------------------------##
####    Ecological Assembly Processes Estimation   ####
##---------------------------------------------------##

# This script was developped by Guillaume Schwob and is based on the scripts developped by Stegen et al. 2013 (https://github.com/stegen/Stegen_etal_ISME_2013) and Richter-Heitmann et al. 2020 (https://github.com/FranzKrah/raup_crick)
# QPE Analysis

# Load packages
library(phyloseq)
library(metagMisc)
library(microbiome)
library(ape)
library(graph4lg)
library(MicEco)
library(reshape2)
library(dplyr)


##### Functions to load ####

## This code is part of the study:
## Stochastic dispersal rather than deterministic selection explains the spatio-temporal distribution of soil bacteria in a temperate grassland
## doi: 10.3389/fmicb.2020.01391
## Franz-Sebastian Krah
## 02 - 27 - 2019

#' @title raup_crick_abu_par
#' @param com community matrix (spXsite)
#' @param reps number of bootstraps
#' @param ncore number of cores (serial: ncore = 1; parallel > 1)
#' @param classic_metric standardizes the metric to range from -1 to 1
#' @param split_ties adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended
#' @details Parallelized version of the Raup-Crick algorithm for "abundance" data (Stegen et al. 2013).
#' Previous code loops over each pairwise community combination;
#' here we randomize the full community matrix and compute Bray-Curtis for the 
#' full matrix and then conduct subsequent Raup-Crick calculations as in Stegen.
#' This increaes computational speed. Further here implemented as multi-core version.
#' The code is acurate and fast (see paper, supplement section D)
#' @author Franz-Sebastian Krah

raup_crick_abu_par <- function(com, reps, ncore, classic_metric=FALSE, split_ties=TRUE){
  
  require("parallel")
  require("doSNOW")
  
  pb <- txtProgressBar(max =reps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  bray.rand <- foreach(randomize = 1:reps, 
                       .options.snow = opts,
                       .packages = c("vegan", "picante")) %dopar% {
                         
                         
                         null.dist <- com*0
                         
                         for(i in 1:nrow(com)){
                           
                           com.pa <- (com>0)*1
                           gamma<-ncol(com)
                           occur<-apply(com>0, MARGIN=2, FUN=sum)
                           abundance<-apply(com, MARGIN=2, FUN=sum)
                           com1 <- rep(0,gamma)
                           
                           com1[sample(1:gamma, sum(com.pa[i,]), replace=FALSE, prob=occur)]<-1
                           com1.samp.sp = sample(which(com1>0), (sum(com[i,])-sum(com1)),
                                                 replace=TRUE,prob=abundance[which(com1>0)]);
                           com1.samp.sp = cbind(com1.samp.sp,1)
                           com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum))
                           colnames(com1.sp.counts) = 'counts'
                           com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
                           com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
                           x <- com1
                           null.dist[i,] <- x
                           rm('com1.samp.sp','com1.sp.counts')
                         }
                         as.matrix(vegdist(null.dist, "bray"))
                       }
  stopCluster(cl)
  
  ## Calculate beta-diversity for obs metacommunity
  bray.obs <- as.matrix(vegdist(com, "bray"))
  
  ##how many null observations is the observed value tied with?
  null_bray_curtis <- bray.rand
  num_exact_matching_in_null <- lapply(null_bray_curtis, function(x) x==bray.obs)
  num_exact_matching_in_null <- apply(simplify2array(num_exact_matching_in_null), 1:2, sum)
  
  ##how many null values are smaller than the observed *dissimilarity*?
  num_less_than_in_null <- lapply(null_bray_curtis, function(x) (x<bray.obs)*1)
  num_less_than_in_null <- apply(simplify2array(num_less_than_in_null), 1:2, sum)
  
  
  rc = (num_less_than_in_null)/reps; # rc;
  
  if(split_ties){
    
    rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
  };
  
  
  if(!classic_metric){
    
    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
    
    rc = (rc-.5)*2
  };
  
  return(rc)
  
}


# Set working directory
setwd("~/GitHub/Habitat_specificity_bacterial_biogeography_SO")

#Load phyloseq
load("physeq_habitat_specificity.RData") 

#### Sediment ####
###### Outputs acquisition ####
#Abundances
sed = subset_samples(physeq_habitat_specificity, Type == "Sediment")
sed.dat = metagMisc::phyloseq_to_df(sed, addtax = F, addtot = F, addmaxrank = F, sorting = "abundance")
rownames(sed.dat) <- sed.dat[,1] #Assign OTU name to the rownames
sed.dat[,1] <- NULL #remove first column
sed.dat=as.data.frame(sed.dat) #transform data frame
sed.dat.t=t(sed.dat) #transpose data frame row=samples, columns=OTUs
meta.sed = microbiome::meta(sed) #Metadata
c<-as.data.frame(meta.sed$Type) #vector with type of habitat
row.names(c) = meta.sed[,1] #assign samples names as rownames 
names(c)[names(c) == 'meta.sed$Type'] <- 'Type'
mer.sed<-merge(c, sed.dat.t, by=0, all=TRUE) #merge vector habitat with abundance data
row.names(mer.sed) = mer.sed$Row.names
com.sed <- mer.sed #community data
env.sed <- com.sed[, 1:2] #information of samples
com.sed <- com.sed[, -c(1:2)] #remove 2 first columns
rowSums(com.sed>0) #total abundances per sample
com.sed <- com.sed[rowSums(com.sed>0) > 1, ] #keep samples with at least one abundances
colSums(com.sed>0) #total abundances per OTUs
com.sed <- com.sed[ , colSums(com.sed>0) >= 1, ] #keep OTUs with at least one abundances
dim(com.sed)
#28 1151
order.sed = as.vector(colnames(com.sed)) #To order the phylogenetic matrix later

# for phylogenetic beta turnover
tree.sed = physeq_habitat_specificity@phy_tree #Extract the tree
tree.sed$tip.label
pruned.phy.sed <- ape::drop.tip(tree.sed, tree.sed$tip.label[!tree.sed$tip.label %in% colnames(com.sed)]) #Subset the tree with the OTUs from sediment
pruned.phy.sed$tip.label
phy.dist.sed <- cophenetic(pruned.phy.sed) #extract cophenetic distances between OTUs : measure of the phylogenetic distance
phy.dist.sed = graph4lg::reorder_mat(phy.dist.sed, order.sed) #To order row and columns of the distance matrix with the same order as OTU table
# Check the order of the distance matrices.
head(rownames(phy.dist.sed))
head(colnames(phy.dist.sed))
head(colnames(com.sed))

##### RUN once and save the outputs #####
# Run betaMNTD/=-betaNTI
nti.sed <- MicEco::ses.comdistnt(com.sed, phy.dist.sed, null.model = "taxa.labels",
                                 runs = 999, iterations = 999, cores = 14)
saveRDS(nti.sed,"nti_abatus2024_all.sediment.RDS")

# Raup-Crick
rc.sed <- raup_crick_abu_par(com.sed, reps = 999, ncore = 14)
rc2.sed = data.frame(as.matrix(rc.sed))
colnames(rc2.sed)=colnames(nti.sed$comdistnt.obs.z)

saveRDS(rc2.sed,"rc2_abatus2024_all.sediment.RDS")

##### Outputs processing #####
#Import RDS
nti_habitat_spec_sed= readRDS("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/QPE Abatus_new/nti_abatus2024_all.sediment.RDS")
rc2_habitat_spec_sed = readRDS("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/QPE Abatus_new/rc2_abatus2024_all.sediment.RDS")

nti.sed<-nti_habitat_spec_sed$comdistnt.obs.z %>% melt() #Extract results and transform matrix, one value for each samples comparison
dim(nti.sed)
rc.sed<-as.matrix(rc2_habitat_spec_sed) %>% melt() #Extract results and transform matrix, one value for each samples comparison
dim(rc.sed)

#Make a matrix with both results and ad column with the corresponding ecological process based on the BNTI and RC values
matrix_bnti_rc_sed <-nti.sed %>% mutate(RC=rc.sed$value) %>% rename(beta_NTI = value)
matrix_bnti_rc_sed$result = NA 
matrix_bnti_rc_sed$beta_NTI = as.numeric(matrix_bnti_rc_sed$beta_NTI) #transform to numeric
matrix_bnti_rc_sed$RC = as.numeric(matrix_bnti_rc_sed$RC) #transform to numeric
matrix_bnti_rc_sed = matrix_bnti_rc_sed %>% 
  mutate(result = ifelse(beta_NTI >= 2, "Variable selection", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC >= 0.95, "Dispersal limitation", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC < 0.95, "Ecological drift", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC <= -0.95, "Homogenizing dispersal", result)) %>%
  mutate(result = ifelse(beta_NTI <= -2, "Homogeneous selection", result))

#absence of homogenizing selection no values of bnti < -2
levels(factor(matrix_bnti_rc_sed$result))
summary(factor(matrix_bnti_rc_sed$result))
View(matrix_bnti_rc_sed)

meta_aba<-physeq_habitat_specificity@sam_data

# Adding the site information to the melted distance dataframe
matrix_bnti_rc_sed$Site1 <- NA
matching_indices <- match(matrix_bnti_rc_sed$Var1, meta_aba$Group)
matrix_bnti_rc_sed$Site1[!is.na(matching_indices)] <- meta_aba$Site[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_sed$Site2 <- NA
matching_indices <- match(matrix_bnti_rc_sed$Var2, meta_aba$Group)
matrix_bnti_rc_sed$Site2[!is.na(matching_indices)] <- meta_aba$Site[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_sed$Site1<- as.character(matrix_bnti_rc_sed$Site1)
matrix_bnti_rc_sed$Site2<- as.character(matrix_bnti_rc_sed$Site2)

# Adding the type information to the melted dataframe
matrix_bnti_rc_sed$Type1 <- NA
matching_indices <- match(matrix_bnti_rc_sed$Var1, meta_aba$Group)
matrix_bnti_rc_sed$Type1[!is.na(matching_indices)] <- meta_aba$Type[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_sed$Type2 <- NA
matching_indices <- match(matrix_bnti_rc_sed$Var2, meta_aba$Group)
matrix_bnti_rc_sed$Type2[!is.na(matching_indices)] <- meta_aba$Type[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_sed$Type1<- as.character(matrix_bnti_rc_sed$Type1)
matrix_bnti_rc_sed$Type2<- as.character(matrix_bnti_rc_sed$Type2)

# Sorting the matrix
sorted_matrix_bnti_rc_sed <- matrix_bnti_rc_sed %>%
  arrange(Site1, Site2, Type1, Type2)

View(matrix_bnti_rc_sed)

######  Making a data frame with the means per comparisons within and among sites  ####
#Remove NA
matrix_sed = na.omit(matrix_bnti_rc_sed)

matrix_sed$Dataset=NA
matrix_sed$Type=NA

matrix_sed_2 <- matrix_sed %>%
  mutate(Type = ifelse(Type1 == "Sediment" & Type2 == "Sediment", "Sediment", Type)) %>%
  mutate(Dataset = ifelse(Site1 == Site2 &
                            Type1 == "Sediment" & Type2 == "Sediment", "intra_site", Dataset)) %>%
  mutate(Dataset = ifelse(Site1 != Site2 &
                            Type1 == "Sediment" & Type2 == "Sediment", "inter_site", Dataset)) %>%
  na.omit()
nrow(matrix_sed_2)
# 756 rows
counts <- matrix_sed_2 %>%
  mutate(Conca = paste( Site1, Site2, sep = "_")) %>%
  select( -Site1, -Site2)%>%
  group_by(Type,Dataset,Conca,result) %>%
  summarise(Count = n(),.groups = "drop")

# Create a reference dataframe with all possible factor levels
all_levels <- expand.grid(Conca=levels(factor(counts$Conca)),
                          result = levels(factor(c("Dispersal limitation","Homogeneous selection","Ecological drift",
                                                   "Variable selection","Homogenizing dispersal"))),
                          Type = levels(factor(counts$Type)),
                          Dataset=levels(factor(counts$Dataset)))
# Merge the counts with all possible levels and fill missing counts with 0
result_df_sed <- all_levels %>%
  left_join(counts, by = c("Type","Dataset","Conca","result")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count))
# View(result_df_sed)
sum(result_df_sed$Count)

# Transform in %
result_df_sed2 <- result_df_sed %>%
  group_by(Type, Dataset, Conca, result) %>%
  summarise(count.tot=sum(Count)) %>%
  ungroup()

sum(result_df_sed2$count.tot)
# 119000 rows
View(result_df_sed2)

# Calculate the percentage within each condition
result_df_sed3 <- result_df_sed2 %>%
  group_by(Type, Dataset, Conca) %>%
  mutate(percent = (count.tot / sum(count.tot)) * 100) %>%
  ungroup() %>%
  na.omit()

sum(result_df_sed3$count.tot)
View(result_df_sed3)
# File ok !
result_df3_sed = result_df_sed3

# Calculate the mean and the se
mean_result_df3_sed <- result_df3_sed %>%
  group_by(Type,Dataset,result) %>%
  summarise(mean_percent = mean(percent),
            se = sd(percent) / sqrt(n()))

View(mean_result_df3_sed)

sed_df_inter_intra = mean_result_df3_sed
View(sed_df_inter_intra)
write.table(sed_df_inter_intra , "sed_df_inter_intra.txt",sep = "\t",row.names = TRUE, col.names = NA)
write.table(result_df3_sed , "sed_df_inter_intra_counts.txt",sep = "\t",row.names = TRUE, col.names = NA)




#### Gut content ####
###### Outputs acquisition ####
#Abundances
cont = subset_samples(physeq_habitat_specificity, Type == "Gut content")
cont.dat = metagMisc::phyloseq_to_df(cont, addtax = F, addtot = F, addmaxrank = F, sorting = "abundance")
rownames(cont.dat) <- cont.dat[,1] #Assign OTU name to the rownames
cont.dat[,1] <- NULL #remove first column
cont.dat=as.data.frame(cont.dat) #transform data frame
cont.dat.t=t(cont.dat) #transpose data frame row=samples, columns=OTUs
meta.cont = microbiome::meta(cont) #Metadata
c<-as.data.frame(meta.cont$Type) #vector with type of habitat
row.names(c) = meta.cont[,1] #assign samples names as rownames 
names(c)[names(c) == 'meta.cont$Type'] <- 'Type'
mer.cont<-merge(c, cont.dat.t, by=0, all=TRUE) #merge vector habitat with abundance data
row.names(mer.cont) = mer.cont$Row.names
com.cont <- mer.cont #community data
env.cont <- com.cont[, 1:2] #information of samples
com.cont <- com.cont[, -c(1:2)] #remove 2 first columns
rowSums(com.cont>0) #total abundances per sample
com.cont <- com.cont[rowSums(com.cont>0) > 1, ] #keep samples with at least one abundances
colSums(com.cont>0) #total abundances per OTUs
com.cont <- com.cont[ , colSums(com.cont>0) >= 1, ] #keep OTUs with at least one abundances
dim(com.cont)
#110 1222
order.cont = as.vector(colnames(com.cont)) #To order the phylogenetic matrix later

# for phylogenetic beta turnover
tree.cont = physeq_habitat_specificity@phy_tree #Extract the tree
tree.cont$tip.label
pruned.phy.cont <- ape::drop.tip(tree.cont, tree.cont$tip.label[!tree.cont$tip.label %in% colnames(com.cont)]) #Subset the tree with the OTUs from contiment
pruned.phy.cont$tip.label
phy.dist.cont <- cophenetic(pruned.phy.cont) #extract cophenetic distances between OTUs : measure of the phylogenetic distance
phy.dist.cont = graph4lg::reorder_mat(phy.dist.cont, order.cont) #To order row and columns of the distance matrix with the same order as OTU table
# Check the order of the distance matrices.
head(rownames(phy.dist.cont))
head(colnames(phy.dist.cont))
head(colnames(com.cont))

##### RUN once and save the outputs #####
# Run betaMNTD/=-betaNTI
nti.cont <- MicEco::ses.comdistnt(com.cont, phy.dist.cont, null.model = "taxa.labels",
                                 runs = 999, iterations = 999, cores = 14)
saveRDS(nti.cont,"nti_abatus2024_all.content.RDS")

# Raup-Crick
rc.cont <- raup_crick_abu_par(com.cont, reps = 999, ncore = 14)
rc2.cont = data.frame(as.matrix(rc.cont))
colnames(rc2.cont)=colnames(nti.cont$comdistnt.obs.z)

saveRDS(rc2.cont,"rc2_abatus2024_all.content.RDS")

##### Outputs processing #####
#Import RDS
nti_habitat_spec_cont= readRDS("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/QPE Abatus_new/nti_abatus2024_all.content.RDS")
rc2_habitat_spec_cont = readRDS("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/QPE Abatus_new/rc2_abatus2024_all.content.RDS")

nti.cont<-nti_habitat_spec_cont$comdistnt.obs.z %>% melt() #Extract results and transform matrix, one value for each samples comparison
dim(nti.cont)
rc.cont<-as.matrix(rc2_habitat_spec_cont) %>% melt() #Extract results and transform matrix, one value for each samples comparison
dim(rc.cont)

#Make a matrix with both results and ad column with the corresponding ecological process bacont on the BNTI and RC values
matrix_bnti_rc_cont <-nti.cont %>% mutate(RC=rc.cont$value) %>% rename(beta_NTI = value)
matrix_bnti_rc_cont$result = NA 
matrix_bnti_rc_cont$beta_NTI = as.numeric(matrix_bnti_rc_cont$beta_NTI) #transform to numeric
matrix_bnti_rc_cont$RC = as.numeric(matrix_bnti_rc_cont$RC) #transform to numeric
matrix_bnti_rc_cont = matrix_bnti_rc_cont %>% 
  mutate(result = ifelse(beta_NTI >= 2, "Variable selection", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC >= 0.95, "Dispersal limitation", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC < 0.95, "Ecological drift", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC <= -0.95, "Homogenizing dispersal", result)) %>%
  mutate(result = ifelse(beta_NTI <= -2, "Homogeneous selection", result))

#absence of homogenizing selection no values of bnti < -2
levels(factor(matrix_bnti_rc_cont$result))
summary(factor(matrix_bnti_rc_cont$result))
View(matrix_bnti_rc_cont)

meta_aba<-physeq_habitat_specificity@sam_data

# Adding the site information to the melted distance dataframe
matrix_bnti_rc_cont$Site1 <- NA
matching_indices <- match(matrix_bnti_rc_cont$Var1, meta_aba$Group)
matrix_bnti_rc_cont$Site1[!is.na(matching_indices)] <- meta_aba$Site[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_cont$Site2 <- NA
matching_indices <- match(matrix_bnti_rc_cont$Var2, meta_aba$Group)
matrix_bnti_rc_cont$Site2[!is.na(matching_indices)] <- meta_aba$Site[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_cont$Site1<- as.character(matrix_bnti_rc_cont$Site1)
matrix_bnti_rc_cont$Site2<- as.character(matrix_bnti_rc_cont$Site2)

# Adding the type information to the melted dataframe
matrix_bnti_rc_cont$Type1 <- NA
matching_indices <- match(matrix_bnti_rc_cont$Var1, meta_aba$Group)
matrix_bnti_rc_cont$Type1[!is.na(matching_indices)] <- meta_aba$Type[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_cont$Type2 <- NA
matching_indices <- match(matrix_bnti_rc_cont$Var2, meta_aba$Group)
matrix_bnti_rc_cont$Type2[!is.na(matching_indices)] <- meta_aba$Type[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_cont$Type1<- as.character(matrix_bnti_rc_cont$Type1)
matrix_bnti_rc_cont$Type2<- as.character(matrix_bnti_rc_cont$Type2)

# Sorting the matrix
sorted_matrix_bnti_rc_cont <- matrix_bnti_rc_cont %>%
  arrange(Site1, Site2, Type1, Type2)

View(matrix_bnti_rc_cont)

######  Making a data frame with the means per comparisons within and among sites  ####
#Remove NA
matrix_cont = na.omit(matrix_bnti_rc_cont)

matrix_cont$Dataset=NA
matrix_cont$Type=NA

matrix_cont_2 <- matrix_cont %>%
  mutate(Type = ifelse(Type1 == "Gut content" & Type2 == "Gut content", "Gut content", Type)) %>%
  mutate(Dataset = ifelse(Site1 == Site2 &
                            Type1 == "Gut content" & Type2 == "Gut content", "intra_site", Dataset)) %>%
  mutate(Dataset = ifelse(Site1 != Site2 &
                            Type1 == "Gut content" & Type2 == "Gut content", "inter_site", Dataset)) %>%
  na.omit()
nrow(matrix_cont_2)
# 756 rows
counts <- matrix_cont_2 %>%
  mutate(Conca = paste( Site1, Site2, sep = "_")) %>%
  select( -Site1, -Site2)%>%
  group_by(Type,Dataset,Conca,result) %>%
  summarise(Count = n(),.groups = "drop")

# Create a reference dataframe with all possible factor levels
all_levels <- expand.grid(Conca=levels(factor(counts$Conca)),
                          result = levels(factor(c("Dispersal limitation","Homogeneous selection","Ecological drift",
                                                   "Variable selection","Homogenizing dispersal"))),
                          Type = levels(factor(counts$Type)),
                          Dataset=levels(factor(counts$Dataset)))
# Merge the counts with all possible levels and fill missing counts with 0
result_df_cont <- all_levels %>%
  left_join(counts, by = c("Type","Dataset","Conca","result")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count))
# View(result_df_cont)
sum(result_df_cont$Count)

# Transform in %
result_df_cont2 <- result_df_cont %>%
  group_by(Type, Dataset, Conca, result) %>%
  summarise(count.tot=sum(Count)) %>%
  ungroup()

sum(result_df_cont2$count.tot)
# 119000 rows
View(result_df_cont2)

# Calculate the percentage within each condition
result_df_cont3 <- result_df_cont2 %>%
  group_by(Type, Dataset, Conca) %>%
  mutate(percent = (count.tot / sum(count.tot)) * 100) %>%
  ungroup() %>%
  na.omit()

sum(result_df_cont3$count.tot)
View(result_df_cont3)
# File ok !
result_df3_cont = result_df_cont3

# Calculate the mean and the se
mean_result_df3_cont <- result_df3_cont %>%
  group_by(Type,Dataset,result) %>%
  summarise(mean_percent = mean(percent),
            se = sd(percent) / sqrt(n()))

View(mean_result_df3_cont)

cont_df_inter_intra = mean_result_df3_cont
View(cont_df_inter_intra)
write.table(cont_df_inter_intra , "cont_df_inter_intra.txt",sep = "\t",row.names = TRUE, col.names = NA)
write.table(result_df3_cont , "cont_df_inter_intra_counts.txt",sep = "\t",row.names = TRUE, col.names = NA)




#### Gut tissue ####
###### Outputs acquisition ####
#Abundances
tiss = subset_samples(physeq_habitat_specificity, Type == "Gut tissue")
tiss.dat = metagMisc::phyloseq_to_df(tiss, addtax = F, addtot = F, addmaxrank = F, sorting = "abundance")
rownames(tiss.dat) <- tiss.dat[,1] #Assign OTU name to the rownames
tiss.dat[,1] <- NULL #remove first column
tiss.dat=as.data.frame(tiss.dat) #transform data frame
tiss.dat.t=t(tiss.dat) #transpose data frame row=samples, columns=OTUs
meta.tiss = microbiome::meta(tiss) #Metadata
c<-as.data.frame(meta.tiss$Type) #vector with type of habitat
row.names(c) = meta.tiss[,1] #assign samples names as rownames 
names(c)[names(c) == 'meta.tiss$Type'] <- 'Type'
mer.tiss<-merge(c, tiss.dat.t, by=0, all=TRUE) #merge vector habitat with abundance data
row.names(mer.tiss) = mer.tiss$Row.names
com.tiss <- mer.tiss #community data
env.tiss <- com.tiss[, 1:2] #information of samples
com.tiss <- com.tiss[, -c(1:2)] #remove 2 first columns
rowSums(com.tiss>0) #total abundances per sample
com.tiss <- com.tiss[rowSums(com.tiss>0) > 1, ] #keep samples with at least one abundances
colSums(com.tiss>0) #total abundances per OTUs
com.tiss <- com.tiss[ , colSums(com.tiss>0) >= 1, ] #keep OTUs with at least one abundances
dim(com.tiss)
#59 1172
order.tiss = as.vector(colnames(com.tiss)) #To order the phylogenetic matrix later

# for phylogenetic beta turnover
tree.tiss = physeq_habitat_specificity@phy_tree #Extract the tree
tree.tiss$tip.label
pruned.phy.tiss <- ape::drop.tip(tree.tiss, tree.tiss$tip.label[!tree.tiss$tip.label %in% colnames(com.tiss)]) #Subset the tree with the OTUs from tissiment
pruned.phy.tiss$tip.label
phy.dist.tiss <- cophenetic(pruned.phy.tiss) #extract cophenetic distances between OTUs : measure of the phylogenetic distance
phy.dist.tiss = graph4lg::reorder_mat(phy.dist.tiss, order.tiss) #To order row and columns of the distance matrix with the same order as OTU table
# Check the order of the distance matrices.
head(rownames(phy.dist.tiss))
head(colnames(phy.dist.tiss))
head(colnames(com.tiss))

##### RUN once and save the outputs #####
# Run betaMNTD/=-betaNTI
nti.tiss <- MicEco::ses.comdistnt(com.tiss, phy.dist.tiss, null.model = "taxa.labels",
                                  runs = 999, iterations = 999, cores = 14)
saveRDS(nti.tiss,"nti_abatus2024_all.tissue.RDS")

# Raup-Crick
rc.tiss <- raup_crick_abu_par(com.tiss, reps = 999, ncore = 14)
rc2.tiss = data.frame(as.matrix(rc.tiss))
colnames(rc2.tiss)=colnames(nti.tiss$comdistnt.obs.z)

saveRDS(rc2.tiss,"rc2_abatus2024_all.tissue.RDS")

##### Outputs processing #####
#Import RDS
nti_habitat_spec_tiss= readRDS("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/QPE Abatus_new/nti_abatus2024_all.tissue.RDS")
rc2_habitat_spec_tiss = readRDS("~/GitHub/Habitat_specificity_bacterial_biogeography_SO/QPE Abatus_new/rc2_abatus2024_all.tissue.RDS")

nti.tiss<-nti_habitat_spec_tiss$comdistnt.obs.z %>% melt() #Extract results and transform matrix, one value for each samples comparison
dim(nti.tiss)
rc.tiss<-as.matrix(rc2_habitat_spec_tiss) %>% melt() #Extract results and transform matrix, one value for each samples comparison
dim(rc.tiss)

#Make a matrix with both results and ad column with the corresponding ecological process batiss on the BNTI and RC values
matrix_bnti_rc_tiss <-nti.tiss %>% mutate(RC=rc.tiss$value) %>% rename(beta_NTI = value)
matrix_bnti_rc_tiss$result = NA 
matrix_bnti_rc_tiss$beta_NTI = as.numeric(matrix_bnti_rc_tiss$beta_NTI) #transform to numeric
matrix_bnti_rc_tiss$RC = as.numeric(matrix_bnti_rc_tiss$RC) #transform to numeric
matrix_bnti_rc_tiss = matrix_bnti_rc_tiss %>% 
  mutate(result = ifelse(beta_NTI >= 2, "Variable selection", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC >= 0.95, "Dispersal limitation", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC < 0.95, "Ecological drift", result)) %>%
  mutate(result = ifelse(beta_NTI < 2 & RC <= -0.95, "Homogenizing dispersal", result)) %>%
  mutate(result = ifelse(beta_NTI <= -2, "Homogeneous selection", result))

#absence of homogenizing selection no values of bnti < -2
levels(factor(matrix_bnti_rc_tiss$result))
summary(factor(matrix_bnti_rc_tiss$result))
View(matrix_bnti_rc_tiss)

meta_aba<-physeq_habitat_specificity@sam_data

# Adding the site information to the melted distance dataframe
matrix_bnti_rc_tiss$Site1 <- NA
matching_indices <- match(matrix_bnti_rc_tiss$Var1, meta_aba$Group)
matrix_bnti_rc_tiss$Site1[!is.na(matching_indices)] <- meta_aba$Site[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_tiss$Site2 <- NA
matching_indices <- match(matrix_bnti_rc_tiss$Var2, meta_aba$Group)
matrix_bnti_rc_tiss$Site2[!is.na(matching_indices)] <- meta_aba$Site[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_tiss$Site1<- as.character(matrix_bnti_rc_tiss$Site1)
matrix_bnti_rc_tiss$Site2<- as.character(matrix_bnti_rc_tiss$Site2)

# Adding the type information to the melted dataframe
matrix_bnti_rc_tiss$Type1 <- NA
matching_indices <- match(matrix_bnti_rc_tiss$Var1, meta_aba$Group)
matrix_bnti_rc_tiss$Type1[!is.na(matching_indices)] <- meta_aba$Type[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_tiss$Type2 <- NA
matching_indices <- match(matrix_bnti_rc_tiss$Var2, meta_aba$Group)
matrix_bnti_rc_tiss$Type2[!is.na(matching_indices)] <- meta_aba$Type[matching_indices[!is.na(matching_indices)]]

matrix_bnti_rc_tiss$Type1<- as.character(matrix_bnti_rc_tiss$Type1)
matrix_bnti_rc_tiss$Type2<- as.character(matrix_bnti_rc_tiss$Type2)

# Sorting the matrix
sorted_matrix_bnti_rc_tiss <- matrix_bnti_rc_tiss %>%
  arrange(Site1, Site2, Type1, Type2)

View(matrix_bnti_rc_tiss)

######  Making a data frame with the means per comparisons within and among sites  ####
#Remove NA
matrix_tiss = na.omit(matrix_bnti_rc_tiss)

matrix_tiss$Dataset=NA
matrix_tiss$Type=NA

matrix_tiss_2 <- matrix_tiss %>%
  mutate(Type = ifelse(Type1 == "Gut tissue" & Type2 == "Gut tissue", "Gut tissue", Type)) %>%
  mutate(Dataset = ifelse(Site1 == Site2 &
                            Type1 == "Gut tissue" & Type2 == "Gut tissue", "intra_site", Dataset)) %>%
  mutate(Dataset = ifelse(Site1 != Site2 &
                            Type1 == "Gut tissue" & Type2 == "Gut tissue", "inter_site", Dataset)) %>%
  na.omit()
nrow(matrix_tiss_2)
# 3422 rows
counts <- matrix_tiss_2 %>%
  mutate(Conca = paste( Site1, Site2, sep = "_")) %>%
  select( -Site1, -Site2)%>%
  group_by(Type,Dataset,Conca,result) %>%
  summarise(Count = n(),.groups = "drop")

# Create a reference dataframe with all possible factor levels
all_levels <- expand.grid(Conca=levels(factor(counts$Conca)),
                          result = levels(factor(c("Dispersal limitation","Homogeneous selection","Ecological drift",
                                                   "Variable selection","Homogenizing dispersal"))),
                          Type = levels(factor(counts$Type)),
                          Dataset=levels(factor(counts$Dataset)))
# Merge the counts with all possible levels and fill missing counts with 0
result_df_tiss <- all_levels %>%
  left_join(counts, by = c("Type","Dataset","Conca","result")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count))
# View(result_df_tiss)
sum(result_df_tiss$Count)

# Transform in %
result_df_tiss2 <- result_df_tiss %>%
  group_by(Type, Dataset, Conca, result) %>%
  summarise(count.tot=sum(Count)) %>%
  ungroup()

sum(result_df_tiss2$count.tot)
# 3422 rows
View(result_df_tiss2)

# Calculate the percentage within each condition
result_df_tiss3 <- result_df_tiss2 %>%
  group_by(Type, Dataset, Conca) %>%
  mutate(percent = (count.tot / sum(count.tot)) * 100) %>%
  ungroup() %>%
  na.omit()

sum(result_df_tiss3$count.tot)
View(result_df_tiss3)
# File ok !
result_df3_tiss = result_df_tiss3

# Calculate the mean and the se
mean_result_df3_tiss <- result_df3_tiss %>%
  group_by(Type,Dataset,result) %>%
  summarise(mean_percent = mean(percent),
            se = sd(percent) / sqrt(n()))

View(mean_result_df3_tiss)

tiss_df_inter_intra = mean_result_df3_tiss
View(tiss_df_inter_intra)
write.table(tiss_df_inter_intra , "tiss_df_inter_intra.txt",sep = "\t",row.names = TRUE, col.names = NA)
write.table(result_df3_tiss , "tiss_df_inter_intra_counts.txt",sep = "\t",row.names = TRUE, col.names = NA)


#### Creating definitive tables ####
#Means within and among sites
ecological_processes_habitat = rbind(sed_df_inter_intra, cont_df_inter_intra, tiss_df_inter_intra)
write.table(ecological_processes_habitat , "ecological_processes_habitat.txt",sep = "\t",row.names = TRUE, col.names = NA)
#Means by sites comparisons
ecological_processes_habitat_counts = rbind(result_df3_sed, result_df3_cont, result_df3_tiss)
write.table(ecological_processes_habitat_counts , "ecological_processes_habitat_counts.txt",sep = "\t",row.names = TRUE, col.names = NA)
