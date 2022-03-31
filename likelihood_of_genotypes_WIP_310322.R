######################################################
###### Code for specific genotypes (SNP by SNP) ######
######################################################
### want to know what the estimated distribution for each loci is
### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c

library(tidyverse)
library(MASS)

read.table("dmi_m0.2_c0.9.main", sep=",", header=T)->data

###only want the 6th deme again - do this up here because the computer is getting cranky about memory
data[which(data$deme==6), ]->data_deme_6
### currently SNPs are wide data, not long. Should I make long so it's on a SNP level not an individual level?

data_long<-gather(data_deme_6, snp, genotype, l1.1:l10.51, factor_key=TRUE)
#data_long$genotype/2->data_long$genotype_scale
#data_long$genotype_scale<-ifelse(data_long$genotype_scale==0, data_long$genotype_scale+0.0001,ifelse(data_long$genotype_scale==1,data_long$genotype_scale-0.0001, data_long$genotype_scale))

### test data to get fitdistr to work

test_data<-data_long[which(data_long$snp=="l1.1"),]
test_data[which(test_data$gen==10), ]->test_data2

nrow<-5100

likelihoods<-matrix(nrow=nrow, ncol=20)
## just for now

####excellent, this seems to work. Am I sure that the size is right (number of individuals *2?)

for(j in 1:20){ ##change this back to 20
  mean_freq<-with(data_long[data_long$rep!=j,], aggregate(genotype, list(gen,snp), mean))
  names(mean_freq)<-c("gen", "snp", "mean_freq")
  mean_freq$index<-paste(mean_freq$gen, mean_freq$snp)#### there is definitely a cleverer way to do this
  data_long$index<-paste(data_long$gen, data_long$snp)
  data_long_leftout<-data_long[data_long$rep==j, ]
  for(i in 1:length(mean_freq$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i, j]<-sum(dbinom(data_long_leftout[data_long_leftout$index==mean_freq$index[i], ]$genotype, size=300, prob=mean_freq$mean_freq[i]/2, log=TRUE))
   }
}
