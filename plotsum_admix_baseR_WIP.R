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

png(file="genomic_cline_plots.png", width=1200, height=1200, units='px')
#png(file="genomic_cline_plots2.4.png", width=1200, height=1200, units='px')
par(mfrow=c(4,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(21,19,20,24,22,23,
                9,7,8,12,10,11,
                3,1,2,6,4,5,
                15,13,14,18,16,17), 4, 6, byrow=TRUE))
for(i in 1:length(unique(summaries_q$index_nosnp))){
plot(0, type="n", xlab="", ylab="", ylim=c(0,1), xlim=c(1,2.5), cex.axis=1.5)
for(j in 1:20){
  lines(summaries_q[which(summaries_q$index_nosnp == unique(summaries_q$index_nosnp)[i]) & summaries_q$rep==j,6], summaries_q[which(summaries_q$index_nosnp == unique(summaries_q$index_nosnp)[i]) & summaries_q$rep==j,8], type='l', col=colours[j])
  }
}

summaries_q$snp_num<2.5 & 
  
  
plot_sum<-ggplot(summaries_q[which(summaries_q$snp_num<2.5),])+geom_line(aes(snp_num, q, colour=as.factor(rep)))+scale_colour_manual(values=colours)+geom_line(aes(snp_num, mean_q), colour='black')+facet_grid(mech~m+c)+theme_bw()
plot_sum<-plot_sum+ylab("Admixture Proportion")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.title.x=element_blank())+guides(fill="none")+theme(legend.position = "none") 

summaries_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'rep','snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_Q$partial_index<-paste(summaries_Q$m, summaries_Q$c, summaries_Q$mech, summaries_Q$snp_num)
summaries_mean_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_Q$partial_index<-paste(summaries_mean_Q$m, summaries_mean_Q$c, summaries_mean_Q$mech, summaries_mean_Q$snp_num)
summaries_mean_Q$index_nosnp<-paste(summaries_mean_Q$m, summaries_mean_Q$c, summaries_mean_Q$mech)
summaries_mean_Q<-summaries_mean_Q[,c(10, 6)]
names(summaries_mean_Q)<-c("partial_index", "mean_Q")
summaries_Q<-merge(summaries_Q,summaries_mean_Q, by='partial_index')

plot_sum_Q<-ggplot(summaries_Q[which(summaries_Q$snp_num<2.5),])+geom_line(aes(snp_num, Q, colour=as.factor(rep)))+scale_colour_manual(values=colours)+geom_line(aes(snp_num, mean_Q), colour='black')+facet_grid(mech~m+c)+theme_bw()+ theme(strip.background = element_blank(),strip.text.x = element_blank())
plot_sum_Q<-plot_sum_Q+xlab("Chromosomes")+ylab("Intersource Ancestry")+guides(colour=guide_legend(title="Replicate"))
plot_sum/plot_sum_Q

ggsave("plotsum_admix_intersource_10.png")