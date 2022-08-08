######################################################
###### Code for specific genotypes (SNP by SNP) ######
######################################################
### want to know what the estimated distribution for each loci is
### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c
#### for this to work, we need to allocate more than default memory on teton eg "salloc --account=project --time=45:00 --nodes=1 --ntasks-per-node=1 --mem=500G"

library(tidyverse)
library(MASS)
library(data.table)

### as of right now, I can't make this happen on my home computer. Needs to be run on Teton
#datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims" , pattern="*main")
### on teton

#### as of right now, I can't even load this into R?

datafiles<-list.files("/gscratch/buerkle/data/incompatible/runs",  pattern="*main", recursive=TRUE, include.dirs=TRUE)
datafiles_11<-datafiles[c(293:316)]
basenames<-basename(datafiles_11)
m<-str_extract(basenames, "(\\d+\\.*\\d*)")
c<-str_match(basenames,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(basenames, "^([^_]+_){1}([^_])") 

alldata<-list()
setwd("/gscratch/buerkle/data/incompatible/runs")

for(i in 1:length(datafiles_11)){
  alldata[[i]]<-fread(datafiles_11[i], sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}

setwd("/gscratch/emcfarl2/predicting_hybrids")
alldata_df<-do.call(rbind.data.frame, alldata)

####clean up
rm(alldata)
#summary(alldata_df)

###only want the 6th deme again - do this up here because the computer is getting cranky about memory
alldata_df[which(alldata_df$deme==6), ]->data_deme_6
### currently SNPs are wide data, not long. Should I make long so it's on a SNP level not an individual level?
#summary(data_deme_6)


data_long<-gather(data_deme_6, snp, genotype, l1.1:l10.51, factor_key=TRUE)

###label the snps that were under selection
data_long[which(data_long$mech %in% c("path_m", "path_e")),]$snp_select<-ifelse(data_long[which(data_long$mech %in% c("path_m", "path_e")),]$snp %in% c("l1.2","l1.20", "l1.4", "l1.40", "l1.6", "l1.60", "l1.8", "l1.80") , 'sel', 'not_sel')
data_long[which(data_long$mech %in% c("dmi_m", "dmi_e")),]$snp_select<-ifelse(data_long[which(data_long$mech %in% c("dmi_m", "dmi_e")),]$snp %in% c("l1.4", "l1.40", "l1.6", "l1.60"), "sel", "not_sel")
data_long[which(data_long$c==0), ]$snp_select<-'not_sel'

##more clean up
rm(alldata_df, data_deme_6)

###only the loci under selection

data_long[which(data_long$snp_select=='sel'),]->data_long_sel
data_long_sel$index<-paste(data_long_sel$m, data_long_sel$c, data_long_sel$gen, data_long_sel$mech, data_long_sel$snp)
nrow_sel<-length(unique(data_long_sel$index))

Fstats_genotypes_sel<-numeric(length = length(unique(data_long_sel$index)))
pvalues_genotypes_sel<-numeric(length(unique(data_long_sel$index)))


for(i in 1:length(unique(data_long_sel$index))){ 
  Fstats_genotypes_sel[[i]]<-unlist(summary(aov(data_long_sel[which(data_long_sel$index==unique(data_long_sel$index)[i]),]$genotype~as.factor(data_long_sel[which(data_long_sel$index==unique(data_long_sel$index)[i]),]$rep))))[7]
  pvalues_genotypes_sel[[i]]<-unlist(summary(aov(data_long_sel[which(data_long_sel$index==unique(data_long_sel$index)[i]),]$genotype~as.factor(data_long_sel[which(data_long_sel$index==unique(data_long_sel$index)[i]),]$rep))))[9]
}

geno_table_sel<-data.frame(model=as.character(unique(data_long_sel$index)), Fstats=Fstats_genotypes_sel, pvalues=pvalues_genotypes_sel)

save.image(file="Predicting_hybrid_genotypes.RData")



#### all loci
data_long$index<-paste(data_long$m, data_long$c, data_long$gen, data_long$mech, data_long$snp) ### took this step out of the loop because it takes a really long time and doesn't need to be done again and again

nrow<-length(unique(data_long$index))


## just for now

####excellent, this seems to work. Am I sure that the size is right (number of individuals *2?)
Fstats_genotypes<-numeric(length = length(unique(data_long$index)))
pvalues_genotypes<-numeric(length(unique(data_long$index)))


for(i in 1:length(unique(data_long$index))){ 
  Fstats_genotypes[[i]]<-unlist(summary(aov(data_long[which(data_long$index==unique(data_long$index)[i]),]$genotype~as.factor(data_long[which(data_long$index==unique(data_long$index)[i]),]$rep))))[7]
 pvalues_genotypes[[i]]<-unlist(summary(aov(data_long[which(data_long$index==unique(data_long$index)[i]),]$genotype~as.factor(data_long[which(data_long$index==unique(data_long$index)[i]),]$rep))))[9]
}

geno_table<-data.frame(model=as.character(unique(data_long$index)), Fstats=Fstats_genotypes, pvalues=pvalues_genotypes)

###clean up before saving!###
rm(data_long, data_long_leftout, data_long_sel, pvalues_genotypes_sel, pvalues_genotypes, Fstats_genotypes, Fstats_genotypes_sel)

save.image(file="Predicting_hybrid_genotypes.RData")
