#### want to fold so that I use 'minor parent' rather than 'raw q' as suggested by Doro

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

alldata_df[which(alldata_df$deme==6), ]->data_deme_6
data_deme_6[which(data_deme_6$gen==10),]->data_deme_6_gen_10
data_deme_6_gen_10[,c(1:8, 519:521, 12, 18, 114)]->alldata_df_6_10
data_long<-gather(alldata_df_6_10, snp, genotype, l1.4:l3.4, factor_key=TRUE)

source('summarySE.R')

data_long$mech<-relevel(data_long$mech, "path_e")
data_long$mech<-relevel(data_long$mech, "path_m")
data_long$mech<-relevel(data_long$mech, "dmi_e")
data_long$mech<-relevel(data_long$mech, "dmi_m")
data_long$mech<-ifelse(data_long$mech=="dmi_m", "dmi", ifelse(data_long$mech=="path_m", "path",ifelse(data_long$mech=="path_e", 'path_e', "dmi_e")))

data_long$snp_num<-(as.numeric(unlist(str_extract(as.factor(data_long$snp),"[[:digit:]]+\\.*[[:digit:]]*"))))
summary(data_long$snp_num)
data_long$q<-data_long$genotype/2
data_long$Q<-ifelse(data_long$genotype==1, 1, 0) 

data_long_noE<-data_long[which(data_long$mech %in% c("dmi", "path")),]
summaries_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'rep'), na.rm=FALSE, conf.interval=.95)
summaries_q$partial_index<-paste(summaries_q$m, summaries_q$c, summaries_q$mech, summaries_q$rep)

data_long_noE$partial_index<-paste(data_long_noE$m, data_long_noE$c, data_long_noE$mech, data_long_noE$rep)

merge(data_long_noE, summaries_q[,c(10, 6)], by='partial_index')->data_long_merged
data_long_merged$rep<-as.factor(data_long_merged$rep)
data_long_merged$minorancestry<-ifelse(data_long_merged$q.y<=0.5, data_long_merged$q.x, 1-data_long_merged$q.x)
colours<-met.brewer(name='OKeeffe1', n=20, type='continuous') 

### want to do 1) correlations between populations and 2) anovas between populations for each SNP
### the schumer papers do this based on haplotypes, I think, not genotypes. I'm going to set it up for genotypes, but then I can do windows? Except we don't have a linkage map...
l1.4_minorancestry_fstat<-vector(length = length(unique(data_long_merged$scenario)))
l1.4_minorancestry_pvalue<-vector(length = length(unique(data_long_merged$scenario)))
correlation_1.4<-vector(length = length(unique(data_long_merged$scenario)))

l1.10_minorancestry_fstat<-vector(length = length(unique(data_long_merged$scenario)))
l1.10_minorancestry_pvalue<-vector(length = length(unique(data_long_merged$scenario)))
correlation_1.10<-vector(length = length(unique(data_long_merged$scenario)))

l3.4_minorancestry_fstat<-vector(length = length(unique(data_long_merged$scenario)))
l3.4_minorancestry_pvalue<-vector(length = length(unique(data_long_merged$scenario)))
correlation_3.4<-vector(length = length(unique(data_long_merged$scenario)))

for(i in 1:length(unique(data_long_merged$scenario))){
  l1.4_minorancestry_fstat[[i]]<-unlist(summary(aov(minorancestry~rep, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4),])))[7]
  l1.4_minorancestry_pvalue[[i]]<-unlist(summary(aov(minorancestry~rep,data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4),])))[9]
  
  l1.10_minorancestry_fstat[[i]]<-unlist(summary(aov(minorancestry~rep, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10),])))[7]
  l1.10_minorancestry_pvalue[[i]]<-unlist(summary(aov(minorancestry~rep,data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10),])))[9]
  
  l3.4_minorancestry_fstat[[i]]<-unlist(summary(aov(minorancestry~rep, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4),])))[7]
  l3.4_minorancestry_pvalue[[i]]<-unlist(summary(aov(minorancestry~rep,data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4),])))[9]
  
  ###let's also do a correlation between the ancestries for rep 1 and 2 for each scenario for each snp
  correlation_1.4[[i]]<-cor.test(data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4 & data_long_merged$rep==1),]$minorancestry, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4 & data_long_merged$rep==2),]$minorancestry)$estimate
  correlation_1.10[[i]]<-cor.test(data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10 & data_long_merged$rep==1),]$minorancestry, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10 & data_long_merged$rep==2),]$minorancestry)$estimate
  correlation_3.4[[i]]<-cor.test(data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4 & data_long_merged$rep==1),]$minorancestry, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4 & data_long_merged$rep==2),]$minorancestry)$estimate
  

  }

minorancestry.Anova<-cbind(unique(as.character(data_long_merged$scenario)), l1.4_minorancestry_fstat, l1.4_minorancestry_pvalue, correlation_1.4, l1.10_minorancestry_fstat, l1.10_minorancestry_pvalue, correlation_1.10, l3.4_minorancestry_fstat, l3.4_minorancestry_pvalue, correlation_3.4)  
write.table(minorancestry.Anova, file="minorancestry.Anova.csv", col.names=T, row.names=F, quote=F, sep=',') 
