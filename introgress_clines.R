##### genomic clines across replicates

### want to look at anovas for SNPs under selection, SNPs on the same chromosome but not under selection and SNPs on a different chromosome, not under selection
library(tidyverse)
library(MASS)
library(data.table)

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
source("genomic.cline.plot.R")
### want to use introgress::clines.plot, I think ###

### might need to build 'cline.data' first, what does this look like?
###let's build alldata_df for deme 6, gen 10
alldata_df_6_10<-alldata_df[which(alldata_df$deme==6 & alldata_df$gen==10),]

###I'm not sure it makes sense using all 510 snps for this. let's do just the 3?
alldata_df_6_10[,c(1:8, 519:521, 12, 18, 63)]->alldata_df_6_10

alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "path_e")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "path_m")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "dmi_e")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "dmi_m")

### I want to do this separately for the 24 categories we have! ###
alldata_df_6_10$index<-paste(alldata_df_6_10$m, alldata_df_6_10$c, alldata_df_6_10$mech)
alldata_df_6_10$index<-as.factor(alldata_df_6_10$index)

genomic.clines<-list() ##oh no, a list!

png(file="genomic_cline_plots.png", width=1200, height=1200, units='px')
par(mfrow=c(4,6),  oma = c(5, 5, 2, 2), mai=c(0.3, 0.3, 0.1, 0.1), mar=c(5,5,1,1))
for(i in length(unique(alldata_df_6_10$index)):1){
  main_label<-unique(alldata_df_6_10$index)[i]
  introgress.data<-alldata_df_6_10[which(alldata_df_6_10$index==unique(alldata_df_6_10$index)[i]),c(12:14)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_6_10[which(alldata_df_6_10$index==unique(alldata_df_6_10$index)[i]),]$q
  loci.data<-matrix(nrow=length(introgress.data), ncol=3)
  dim(loci.data)<-c(3, 3)
  loci.data[,1]<-colnames(introgress.data)
  loci.data[,2]<-'C'
  loci.data[,3]<-chromosome
  genomic.clines[[i]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=alldata_df_6_10[which(alldata_df_6_10$index==unique(alldata_df_6_10$index)[i]),]$q, loci.data=loci.data)
  genomic.cline.plot(genomic.clines[[i]])
}
mtext('Admixture Proportion', side = 1, outer = TRUE, line = 2, cex=2.5)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2.5)
dev.off()

#### need to use this to get the order I want!!
layout(matrix(c(25,25,25,25,1,2,3,4,5,6,
                25,25,25,25,7,8,9,10,11,12,
                25,25,25,25,13,14,15,16,17,18,
                25,25,25,25,19,20,21,22,23,24), 4, 10, byrow=TRUE))