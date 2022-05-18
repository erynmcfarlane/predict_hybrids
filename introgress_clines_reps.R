##### genomic clines across replicates

### want to look at anovas for SNPs under selection, SNPs on the same chromosome but not under selection and SNPs on a different chromosome, not under selection
library(tidyverse)
library(MASS)
library(data.table)
library(MetBrewer)
library(introgress)

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
#save.image("introgress_working.RData")
#load("introgress_working.RData")

source("genomic.cline.plot.reps.R")
### want to use introgress::clines.plot, I think ###

### might need to build 'cline.data' first, what does this look like?
###let's build alldata_df for deme 6, gen 10
alldata_df_6_10<-alldata_df[which(alldata_df$deme==6 & alldata_df$gen==10),]

###I'm not sure it makes sense using all 510 snps for this. let's do just the 3?
alldata_df_6_10[,c(1:8, 519:521, 12, 18, 63)]->alldata_df_6_10
alldata_df_6_10$mech<-as.factor(ifelse(alldata_df_6_10$mech=="dmi_m", "dmi", ifelse(alldata_df_6_10$mech=="path_m", "path",ifelse(alldata_df_6_10$mech=="path_e", 'path_e', "dmi_e"))))

alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "path_e")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "path")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "dmi_e")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "dmi")

### I want to do this separately for the 24 categories we have! ###
alldata_df_6_10$index<-paste(alldata_df_6_10$m, alldata_df_6_10$c, alldata_df_6_10$mech)
alldata_df_6_10$index_reps<-paste(alldata_df_6_10$m, alldata_df_6_10$c, alldata_df_6_10$mech, alldata_df_6_10$rep)
alldata_df_6_10$index<-as.factor(alldata_df_6_10$index)

alldata_df_6_10_noE<-alldata_df_6_10[which(alldata_df_6_10$mech %in% c("dmi", "path")),]

##### this is really close, just keep on it!
genomic.clines<-list() ##oh no, a list!
colours<-met.brewer(name='OKeeffe1', n=20, type='continuous') 
#png(file="genomic_cline_plots1.4.png", width=1200, height=1200, units='px')
#png(file="genomic_cline_plots1.10.png", width=1200, height=1200, units='px')
png(file="genomic_cline_plots2.4.png", width=1200, height=1200, units='px')
par(mfrow=c(2,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(9,7,8,12,10,11,
                3,1,2,6,4,5), 2, 6, byrow=TRUE))
for(i in 1:length(unique(alldata_df_6_10_noE$index))){
  genomic.clines.reps<-list()
    for(j in 1:20){
  introgress.data<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),c(12:14)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),]$q
  loci.data<-matrix(nrow=length(introgress.data), ncol=3)
  dim(loci.data)<-c(3, 3)
  loci.data[,1]<-colnames(introgress.data)
  loci.data[,2]<-'C'
  loci.data[,3]<-chromosome
  genomic.clines.reps[[j]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=hi.index, loci.data=loci.data)
    }
  plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), cex.axis=1.5)
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
  genomic.cline.plot(genomic.clines.reps)
  if (i %in% c(9,7,8,12,10,11))
    {
    mtext(paste0("m = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$m[i]), line=2, cex=1.5)
    mtext(paste0("c = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$c[i]), line=0.25, cex=1.5)
    }
  
  if (i %in% c(11,5))
  {
    if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
    }
}

mtext('Admixture Proportion', side = 1, outer = TRUE, line = 2, cex=2.5)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2.5)
dev.off()

