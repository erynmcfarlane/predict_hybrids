##### genomic clines across replicates

### want to look at anovas for SNPs under selection, SNPs on the same chromosome but not under selection and SNPs on a different chromosome, not under selection
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

### skip down for multinomial clines (line 209) ####

library(introgress)

### want to use introgress::clines.plot, I think ###

### might need to build 'cline.data' first, what does this look like?
###let's build alldata_df for deme 6, gen 10
alldata_df_6_10<-alldata_df[which(alldata_df$deme==6 & alldata_df$gen==10),]

### I want to do this separately for the 24 categories we have! ###
alldata_df_6_10$index<-paste(alldata_df_6_10$m, alldata_df_6_10$c, alldata_df_6_10$mech)
genomic.clines<-list() ##oh no, a list!

#for(i in 1:length(unique(alldata_df_6_10$index))){
for(i in 1:3){
  introgress.data<-alldata_df_6_10[which(alldata_df_6_10$index==unique(alldata_df_6_10$index)[i]),c(9:518)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_6_10[which(alldata_df_6_10$index==unique(alldata_df_6_10$index)[i]),]$q
  loci.data<-matrix(nrow=length(introgress.data), ncol=3)
  dim(loci.data)<-c(510, 3)
  loci.data[,1]<-colnames(introgress.data)
  loci.data[,2]<-'C'
  loci.data[,3]<-chromosome
  genomic.clines[[i]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=alldata_df_6_10[which(alldata_df_6_10$index==unique(alldata_df_6_10$index)[i]),]$q, loci.data=loci.data)
  clines.plot<-composite.clines(genomic.clines[[i]], pdf=FALSE)
  print(png(clines.plot, file=paste("clines, [i], .png"))
}
clines.plot<-composite.clines(genomic.clines[[i]], pdf=FALSE)
png(clines.plot, file="clines_plot_6_10.png")
