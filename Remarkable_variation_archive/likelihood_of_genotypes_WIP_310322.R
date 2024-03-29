######################################################
###### Code for specific genotypes (SNP by SNP) ######
######################################################
### want to know what the estimated distribution for each loci is
### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c
#### for this to work, we need to allocate more than default memory on teton eg "salloc --account=project --time=45:00 --nodes=1 --ntasks-per-node=1 --mem=500G"

library(tidyverse)
library(MASS)
library(data.table)

### as of right now, I can't make this happen on my home computer. Needs to be run on Teton
#datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims" , pattern="*main")
### on teton

#### as of right now, I can't even load this into R?

datafiles<-list.files("/gscratch/buerkle/data/incompatible/runs/deme11",  pattern="*main")
m<-str_extract(datafiles, "(\\d+\\.*\\d*)")
c<-str_match(datafiles,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(datafiles, "^([^_]+_){1}([^_])") 

alldata<-list()
setwd("/gscratch/buerkle/data/incompatible/runs/deme11")

for(i in 1:length(datafiles)){
  alldata[[i]]<-fread(datafiles[i], sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}

setwd("/gscratch/emcfarl2/predicting_hybrids")
alldata_df<-do.call(rbind.data.frame, alldata)

####clean up
rm(alldata)
#summary(alldata_df)

###only want the 6th deme again - do this up here because the computer is getting cranky about memory
alldata_df[which(alldata_df$deme==6), ]->data_deme_6
### currently SNPs are wide data, not long. Should I make long so it's on a SNP level not an individual level?
#summary(data_deme_6)


data_long<-gather(data_deme_6, snp, genotype, l1.1:l10.51, factor_key=TRUE)

##more clean up
rm(alldata_df, data_deme_6)
#data_long$genotype/2->data_long$genotype_scale
#data_long$genotype_scale<-ifelse(data_long$genotype_scale==0, data_long$genotype_scale+0.0001,ifelse(data_long$genotype_scale==1,data_long$genotype_scale-0.0001, data_long$genotype_scale))

### test data to make sure this all works

#test_data<-data_long[which(data_long$snp=="l1.1"),]
#test_data[which(test_data$gen==10), ]->test_data2

data_long$index<-paste(data_long$m, data_long$c, data_long$gen, data_long$mech, data_long$snp) ### took this step out of the loop because it takes a really long time and doesn't need to be done again and again

nrow<-length(unique(data_long$index))

likelihoods<-matrix(nrow=nrow, ncol=20)
## just for now

####excellent, this seems to work. Am I sure that the size is right (number of individuals *2?)

start_time<-Sys.time()

for(j in 1:20){ ##change this back to 20
  mean_freq<-with(data_long[data_long$rep!=j,], aggregate(genotype, list(m, c, gen, mech, snp), mean))
  names(mean_freq)<-c("m", "c", "gen", "mech", "snp", "mean_freq")
  mean_freq$index<-paste(mean_freq$m, mean_freq$c, mean_freq$gen, mean_freq$mech, mean_freq$snp)#### there is definitely a cleverer way to do this
  data_long_leftout<-data_long[data_long$rep==j, ]
  for(i in 1:length(mean_freq$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i, j]<-sum(dbinom(data_long_leftout[data_long_leftout$index==mean_freq$index[i], ]$genotype, size=300, prob=mean_freq$mean_freq[i]/2, log=TRUE))
  }
  print(j)
}

end_time<-Sys.time()
end_time-start_time

###clean up before saving!###
rm(data_long, data_long_leftout)

save.image(file="Predicting_hybrid_genotypes.RData")
