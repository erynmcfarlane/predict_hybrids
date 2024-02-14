#### Input script for SNP data analyses 
##### genomic clines across replicates

### want to look at anovas for SNPs under selection, SNPs on the same chromosome but not under selection and SNPs on a different chromosome, not under selection
library(tidyverse)
library(MASS)
library(data.table)
library(MetBrewer)
library(introgress)

#datafiles<-list.files("/gscratch/buerkle/data/incompatible/runs",  pattern="*main", recursive=TRUE, include.dirs=TRUE)
###folder below is hardcoded for each 2:1, 5:1, 10:1 and 50:1
datafiles<-list.files("/gscratch/emcfarl2/predicting_hybrids/twoone",  pattern="*main", recursive=TRUE, include.dirs=TRUE)
#datafiles_11<-datafiles[c(293:316)]
basenames<-basename(datafiles)
m<-str_extract(basenames, "(\\d+\\.*\\d*)")
c<-str_match(basenames,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(basenames, "^([^_]+_){1}([^_])") 

alldata<-list()
#setwd("/gscratch/buerkle/data/incompatible/runs")

for(i in 1:length(datafiles)){
  alldata[[i]]<-fread(datafiles[i], sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}

#setwd("/gscratch/emcfarl2/predicting_hybrids")
alldata_df<-do.call(rbind.data.frame, alldata)

####clean up
rm(alldata)

### skip down for multinomial clines (line 209) ####
#save.image("introgress_working.RData")
#load("introgress_working.RData")

alldata_df[which(alldata_df$deme==2), ]->data_deme_2
data_deme_2[which(data_deme_2$gen==10),]->data_deme_2_gen_10



