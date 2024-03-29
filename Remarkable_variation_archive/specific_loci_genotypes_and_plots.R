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
### skip to the bottom from here for doro's plots (line 139) ####

### skip down for multinomial clines (line 209) ####
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



source('summarySE.R')

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

unlist(summary(aov(data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.4'), ]$genotype~data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.4'), ]$rep)))[7]
unlist(summary(aov(data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.10'), ]$genotype~data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.4'), ]$rep)))[7]
unlist(summary(aov(data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l2.4'), ]$genotype~data_long[which(data_long$deme==6 & data_long$gen==100 & data_long$snp=='l1.4'), ]$rep)))[7]




snp_sel_6<-data_long[which(data_long$deme==6 & data_long$snp=='l1.4'), ]
snp_ld_6<-data_long[which(data_long$deme==6 & data_long$snp=='l1.10'), ]
snp_nosel_6<-data_long[which(data_long$deme==6 &  data_long$snp=='l2.4'), ]

anova(lm(genotype~rep+mech+c+m+gen, data=snp_sel_6))
anova(lm(genotype~rep+mech+c+m+gen, data=snp_ld_6))
anova(lm(genotype~rep+mech+c+m+gen, data=snp_nosel_6))

anova<-aov(genotype~snp+as.factor(rep)+snp:as.factor(rep)+gen+snp:gen+mech+c+m, data=data_long[which(data_long$deme==6),])
###think a little bit about what exactly you want in here!
Tukey<-TukeyHSD(anova)
#Analysis of Variance Table
#### results from when I ran this on teton
summary(anova)
#Df  Sum Sq Mean Sq F value  Pr(>F)    
#snp                      2      46    22.8   40.66 < 2e-16 ***
 # as.factor(rep)          19    6056   318.7  568.39 < 2e-16 ***
#  gen                      1     104   104.3  185.94 < 2e-16 ***
 # mech                     3    1169   389.7  694.88 < 2e-16 ***
#  as.factor(c)             2    2071  1035.7 1846.94 < 2e-16 ***
  as.factor(m)             1    1830  1830.4 3264.21 < 2e-16 ***
  snp:as.factor(rep)      38    2414    63.5  113.29 < 2e-16 ***
  snp:gen                  2      17     8.4   15.05 2.9e-07 ***
  Residuals          2159931 1211172     0.6    



### might want to look at specific cases - tukey?
###double check that I only have deme 6!


#####GENERATION 10#####
alldata_df[which(alldata_df$deme==6), ]->data_deme_6
data_deme_6[which(data_deme_6$gen==10),]->data_deme_6_gen_10
data_long<-gather(data_deme_6_gen_10, snp, genotype, l1.1:l10.51, factor_key=TRUE)


#### this puts each of the mechanisms in the appropriate order for the plots###
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
summaries_mean_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_q$partial_index<-paste(summaries_mean_q$m, summaries_mean_q$c, summaries_mean_q$mech, summaries_mean_q$snp_num)
summaries_mean_q<-summaries_mean_q[,c(10, 6)]
names(summaries_mean_q)<-c("partial_index", "mean_q")
summaries_q<-merge(summaries_q,summaries_mean_q, by='partial_index')
plot_sum<-ggplot(summaries_q[which(summaries_q$snp_num<2.5),])+geom_line(aes(snp_num, q, colour=as.factor(rep)))+scale_colour_manual(values=colours)+geom_line(aes(snp_num, mean_q), colour='black')+facet_grid(mech~m+c)+theme_bw()
plot_sum<-plot_sum+ylab("Admixture Proportion")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.title.x=element_blank())+guides(fill="none")+theme(legend.position = "none") 
                                                            
summaries_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'rep','snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_Q$partial_index<-paste(summaries_Q$m, summaries_Q$c, summaries_Q$mech, summaries_Q$snp_num)
summaries_mean_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_Q$partial_index<-paste(summaries_mean_Q$m, summaries_mean_Q$c, summaries_mean_Q$mech, summaries_mean_Q$snp_num)
summaries_mean_Q<-summaries_mean_Q[,c(10, 6)]
names(summaries_mean_Q)<-c("partial_index", "mean_Q")
summaries_Q<-merge(summaries_Q,summaries_mean_Q, by='partial_index')

plot_sum_Q<-ggplot(summaries_Q[which(summaries_Q$snp_num<2.5),])+geom_line(aes(snp_num, Q, colour=as.factor(rep)))+scale_colour_manual(values=colours)+geom_line(aes(snp_num, mean_Q), colour='black')+facet_grid(mech~m+c)+theme_bw()+ theme(strip.text.x = element_blank())
plot_sum_Q<-plot_sum_Q+xlab("Chromosomes")+ylab("Intersource Ancestry")+guides(fill="none")+theme(legend.position = "none") 
plot_sum/plot_sum_Q

ggsave("plotsum_admix_intersource_10.png")
