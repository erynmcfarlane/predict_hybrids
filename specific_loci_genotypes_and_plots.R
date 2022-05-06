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
### skip to the bottom from here for doro's plots
alldata_df[,c(1:8, 519:521, 12, 18, 63)]->alldata_df_3snps

data_long<-gather(alldata_df_3snps, snp, genotype, l1.4:l2.4, factor_key=TRUE)

data_long$index<-paste(data_long$deme, data_long$m, data_long$c, data_long$gen, data_long$mech)

length(unique(data_long[which(data_long$deme==6),]$index))

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

### do this again, but do densities of each loci across 3x24 panel plots

snp_sel_densities<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.4'), ], aes(genotype, colour=as.factor(rep)))+geom_density()+facet_grid(mech~m+c)+theme_bw()
ggsave("snp_sel_densities10.png")
snp_ld_densities<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.10'), ], aes(genotype, colour=as.factor(rep)))+geom_density()+facet_grid(mech~m+c)+theme_bw()
ggsave("snp_ld_densities10.png")
snp_nosel_densities<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l2.4'), ], aes(genotype, colour=as.factor(rep)))+geom_density()+facet_grid(mech~m+c)+theme_bw()
ggsave("snp_nosel_densities10.png")


snp_sel_densities<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.4'), ], aes(genotype, colour=as.factor(rep)))+geom_density()+facet_grid(mech~m+c)+theme_bw()
ggsave("snp_sel_densities100.png")
snp_ld_densities<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.10'), ], aes(genotype, colour=as.factor(rep)))+geom_density()+facet_grid(mech~m+c)+theme_bw()
ggsave("snp_ld_densities100.png")
snp_nosel_densities<-ggplot(data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l2.4'), ], aes(genotype, colour=as.factor(rep)))+geom_density()+facet_grid(mech~m+c)+theme_bw()
ggsave("snp_nosel_densities100.png")


###want to compare F for each of these snps compare within 

unlist(summary(aov(data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.4'), ]$genotype~data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.4'), ]$rep)))[7]
unlist(summary(aov(data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.10'), ]$genotype~data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.4'), ]$rep)))[7]
unlist(summary(aov(data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l2.4'), ]$genotype~data_long[which(data_long$deme==6 & data_long$gen==10 & data_long$snp=='l1.4'), ]$rep)))[7]

snp_sel_6<-data_long[which(data_long$deme==6 & data_long$snp=='l1.4'), ]
snp_ld_6<-data_long[which(data_long$deme==6 & data_long$snp=='l1.10'), ]
snp_nosel_6<-data_long[which(data_long$deme==6 &  data_long$snp=='l2.4'), ]

anova(lm(genotype~rep+mech+c+m+gen, data=snp_sel_6))
anova(lm(genotype~rep+mech+c+m+gen, data=snp_ld_6))
anova(lm(genotype~rep+mech+c+m+gen, data=snp_nosel_6))

anova(lm(genotype~snp+rep+snp:rep+mech+c+m+gen, data=data_long))

#Analysis of Variance Table
#### results from when I ran this on teton
#Response: genotype
#Df   Sum Sq Mean Sq  F value    Pr(>F)    
#snp              2       29   14.57  17.1378 3.607e-08 ***
#rep              1        2    2.19   2.5737    0.1087    
#mech             3     1312  437.43 514.4983 < 2.2e-16 ***
# c                1       79   78.77  92.6429 < 2.2e-16 ***
# m                1       98   97.81 115.0399 < 2.2e-16 ***
#gen              1       29   28.87  33.9605 5.624e-09 ***
#snp:rep          2      111   55.45  65.2163 < 2.2e-16 ***
#Residuals 23759988 20200945    0.85


### might want to look at specific cases - tukey?
###double check that I only have deme 6!
Fstats<-list()
pvalues<-list()
for(i in 1:length(unique(data_long[which(data_long$gen==10 & data_long$deme==6),]$index))){
  Fstats[[i]]<-anova(lm(genotype~snp*rep, data=data_long[data_long$index==unique(data_long$index)[i],]))[3,4]
  pvalues[[i]]<-anova(lm(genotype~snp*rep, data=data_long[data_long$index==unique(data_long$index)[i],]))[3,5]
}




#### want to replicate figure 3, A, top and middle rows from Lindtke and Buerkle 2015
### q= genotype/2
### Q = genotype = 1 = 1, if genotype=0 or 2, = 0

alldata_df[which(alldata_df$deme==6), ]->data_deme_6
data_deme_6[which(data_deme_6$gen==10),]->data_deme_6_gen_10
data_long<-gather(data_deme_6_gen_10, snp, genotype, l1.1:l10.51, factor_key=TRUE)


#### this puts each of the mechanisms in the appropriate order for the plots###
data_long$mech<-relevel(data_long$mech, "path_e")
data_long$mech<-relevel(data_long$mech, "path_m")
data_long$mech<-relevel(data_long$mech, "dmi_e")
data_long$mech<-relevel(data_long$mech, "dmi_m")

data_long$snp_num<-(as.numeric(unlist(str_extract(as.factor(data_long$snp),"[[:digit:]]+\\.*[[:digit:]]*"))))

data_long$q<-data_long$genotype/2
### still didn't work. Don't know why
plot1<-ggplot(data_long, aes(q~snp_num, colour=rep))+geom_line()+facet_grid(mech~m+c)+theme_bw()

ggsave("plot1.png")