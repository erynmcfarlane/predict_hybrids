######################################################
###### Code for specific genotypes (SNP by SNP) ######
######################################################
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

alldata_df[,c(1:8, 519:521, 12, 18, 63)]->alldata_df_3snps

data_long<-gather(alldata_df_3snps, snp, genotype, l1.4:l2.4, factor_key=TRUE)

data_long$index<-paste(data_long$deme, data_long$m, data_long$c, data_long$gen, data_long$mech)

length(unique(data_long[which(data_long$deme==6)],$index))

#### this puts each of the mechanisms in the appropriate order for the plots###
data_long$mech<-relevel(data_long$mech, "path_e")
data_long$mech<-relevel(data_long$mech, "path_m")
data_long$mech<-relevel(data_long$mech, "dmi_e")
data_long$mech<-relevel(data_long$mech, "dmi_m")


### could loop this through all of the demes and all of the generations###
deme6_gen10<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==10), ], aes(as.factor(genotype),q, colour=snp))+geom_boxplot()+coord_flip()+facet_grid(mech~m+c)+theme_bw()

ggsave("deme6_gen10.png")

### do this again for generation 100
deme6_gen100<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==100), ], aes(as.factor(genotype),q, colour=snp))+geom_boxplot()+coord_flip()+facet_grid(mech~m+c)+theme_bw()

ggsave("deme6_gen100.png")



source(summarySE.R)

###this shows each individual genotype, but doesn't show the variation among replicates
##happy function that gives me that gives me means and SEs for each of my groups.
summaries<-summarySE(data_long[which(data_long$deme==6 & data_long$gen==10), ], measurevar='q', groupvars=c("m", "c", "mech", 'rep','snp','genotype'), na.rm=FALSE, conf.interval=.95)
deme6_gen10_means<-ggplot(summaries, aes(as.factor(genotype),q, colour=snp))+geom_boxplot()+coord_flip()+facet_grid(mech~m+c)+theme_bw()
ggsave("deme6_gen10_means.png")

summaries<-summarySE(data_long[which(data_long$deme==6 & data_long$gen==100), ], measurevar='q', groupvars=c("m", "c", "mech", 'rep','snp','genotype'), na.rm=FALSE, conf.interval=.95)
deme6_gen100_means<-ggplot(summaries, aes(as.factor(genotype),q, colour=snp))+geom_boxplot()+coord_flip()+facet_grid(mech~m+c)+theme_bw()
ggsave("deme6_gen100_means.png")



### do this again, but do densities of each loci accross 3x24 panel plots

