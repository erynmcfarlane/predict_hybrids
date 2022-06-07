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
#summary(alldata_df)
### skip to the bottom from here for doro's plots (line 139) ####
alldata_df[which(alldata_df$deme==6), ]->data_deme_6
data_deme_6[which(data_deme_6$gen==10),]->data_deme_6_gen_10
data_long<-gather(data_deme_6_gen_10, snp, genotype, l1.1:l10.51, factor_key=TRUE)

source('summarySE.R')

data_long$mech<-relevel(data_long$mech, "path_e")
data_long$mech<-relevel(data_long$mech, "path_m")
data_long$mech<-relevel(data_long$mech, "dmi_e")
data_long$mech<-relevel(data_long$mech, "dmi_m")
data_long$mech<-ifelse(data_long$mech=="dmi_m", "dmi", ifelse(data_long$mech=="path_m", "path",ifelse(data_long$mech=="path_e", 'path_e', "dmi_e")))

data_long$snp_num<-(as.numeric(unlist(str_extract(as.factor(data_long$snp),"[[:digit:]]+\\.*[[:digit:]]*"))))
summary(data_long$snp_num)
data_long$q<-data_long$genotype/2
data_long$Q<-ifelse(data_long$genotype==1, 1, 0) ### is this right? There are really only two locus-specific options, right?

library(MetBrewer)
library(patchwork)
colours=met.brewer(name="OKeeffe1", n=20, type="continuous")
data_long_noE<-data_long[which(data_long$mech %in% c("dmi", "path")),]
summaries_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'rep','snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_q$partial_index<-paste(summaries_q$m, summaries_q$c, summaries_q$mech, summaries_q$snp_num)
summaries_q$index_nosnp<-paste(summaries_q$m, summaries_q$c, summaries_q$mech)
summaries_mean_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_q$partial_index<-paste(summaries_mean_q$m, summaries_mean_q$c, summaries_mean_q$mech, summaries_mean_q$snp_num)
summaries_mean_q<-summaries_mean_q[,c(10, 6)]
names(summaries_mean_q)<-c("partial_index", "mean_q")
summaries_q<-merge(summaries_q,summaries_mean_q, by='partial_index')

summaries_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'rep','snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_Q$partial_index<-paste(summaries_Q$m, summaries_Q$c, summaries_Q$mech, summaries_Q$snp_num)
summaries_Q$index_nosnp<-paste(summaries_Q$m, summaries_Q$c, summaries_Q$mech)
summaries_mean_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_Q$partial_index<-paste(summaries_mean_Q$m, summaries_mean_Q$c, summaries_mean_Q$mech, summaries_mean_Q$snp_num)
summaries_mean_Q<-summaries_mean_Q[,c(10, 6)]
names(summaries_mean_Q)<-c("partial_index", "mean_Q")
summaries_Q<-merge(summaries_Q,summaries_mean_Q, by='partial_index')

summaries_q<-summaries_q[which(summaries_q$snp_num<2.51),] ### to only use the chromosomes I'm gonna plot
summaries_Q<-summaries_Q[which(summaries_Q$snp_num<2.51),]

###admix
pdf(file="plotsum_admix_10.pdf", width=10, height=5)
par(mfrow=c(2,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(1,3,5,7,9,11,
                2,4,6,8,10,12), 2, 6, byrow=TRUE))
for(i in 1:length(unique(summaries_q$index_nosnp))){
plot(0, type="l", xlab="", ylab="", ylim=c(0,1), xlim=c(1,2.5), cex.axis=2)
summaries_q_temp<-summaries_q[which(summaries_q$index_nosnp == unique(summaries_q$index_nosnp)[i]),]
for(j in 1:20){
  lines(summaries_q_temp[which(summaries_q_temp$rep==j),6] , summaries_q_temp[which(summaries_q_temp$rep==j),8],type='l', col=colours[j])
  lines(summaries_q_temp[which(summaries_q_temp$rep==j),6], summaries_q_temp[which(summaries_q_temp$rep==j),13], type='l', col='black')
  }


if (i %in% c(1,3,5,7,9,11))
{
  mtext(paste0("m = ", summaries_q[which(summaries_q$index==unique(summaries_q$index)[i]),]$m[i]), line=2, cex=1.5)
  mtext(paste0("c = ", summaries_q[which(summaries_q$index==unique(summaries_q$index)[i]),]$c[i]), line=0.25, cex=1.5)
}

if (i %in% c(11,12))
{
  if (summaries_q[which(summaries_q$index==unique(summaries_q$index)[i]),]$mech[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
  else if (summaries_q[which(summaries_q$index==unique(summaries_q$index)[i]),]$mech[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
  }
}
mtext('Admixture Proportion', side = 1, outer = TRUE, line = 2, cex=2)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2)
dev.off()

### Intersource Ancestry ###
pdf(file="plotsum_intersource_10.pdf", width=10, height=5)
par(mfrow=c(2,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(1,3,5,7,9,11,
                2,4,6,8,10,12), 2, 6, byrow=TRUE))

for(i in 1:length(unique(summaries_Q$index_nosnp))){
  plot(0, type="l", xlab="", ylab="", ylim=c(0,1), xlim=c(1,2.5), cex.axis=2)
  summaries_Q_temp<-summaries_Q[which(summaries_Q$index_nosnp == unique(summaries_Q$index_nosnp)[i]),]
  for(j in 1:20){
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6] , summaries_Q_temp[which(summaries_Q_temp$rep==j),8],type='l', col=colours[j])
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6], summaries_Q_temp[which(summaries_Q_temp$rep==j),13], type='l', col='black')
  }
  
  
  if (i %in% c(1,3,5,7,9,11))
  {
   mtext(paste0("m = ", summaries_Q[which(summaries_Q$index==unique(summaries_Q$index)[i]),]$m[i]), line=2, cex=1.5)
    mtext(paste0("c = ", summaries_Q[which(summaries_Q$index==unique(summaries_Q$index)[i]),]$c[i]), line=0.25, cex=1.5)
  }
  
  if (i %in% c(11,12))
  {
    if (summaries_Q[which(summaries_Q$index==unique(summaries_Q$index)[i]),]$mech[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
    else if (summaries_Q[which(summaries_Q$index==unique(summaries_Q$index)[i]),]$mech[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
  }
}
mtext('Intersource Ancestry', side = 1, outer = TRUE, line = 2, cex=2)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2)
dev.off()