###want to extend the njunctions loop to het scores

#### the plan is to see if we can use 19 demes from a simulation to predict the 20th, accross different generations. 
###when working on local computer
setwd("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/deme11")

##########################
#### load in the data#####
##########################
library(stringr)
library(MASS)
library(fitdistrplus)
##if using local computer
datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/deme11" , pattern="*first8cols.txt.gz")
##if using teton
#datafiles<-list.files("/project/evolgen/jjahner/hybrid_sims/" , pattern="*first8cols.txt.gz")
m<-str_extract(datafiles, "(\\d+\\.*\\d*)")
c<-str_match(datafiles,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?


### just for the moment because #4 and #5 aren't perfect
datafiles<-datafiles[c(1:3, 6)] ### take this out when we've fixed the grep

alldata<-list()
for(i in 1:length(datafiles)){
  alldata[[i]]<-read.table(gzfile(datafiles[i]), sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]])))
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
}

alldata_df<-do.call(rbind.data.frame, alldata)

summary(alldata_df)

### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)

### only going to use middle, hybridizing demes 5 and 6
## two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))

### this is because the beta distribution is 0 to 1, excluding 0 and 1. I'm not entirely sure how much this matters, but I think it's what's needed to make fdistr beta to work
alldata_df_6$het<-ifelse(alldata_df_6$het ==0, alldata_df_6$het+0.001, ifelse(alldata_df_6$het == 1, alldata_df_6$het-0.001, alldata_df_6$het))

get_dist<-function(x){
  estimates<-fitdistr(x, 'beta', start=list(shape1=10,shape2=10), lower=c(0,0)) ### not at all sure how to decide on the starting shapes!
  estimates$estimate
  ### this completely ignores the error around the alpha and beta distributions. Not sure about this at all. 
}

for(j in 1:20){ ##change this back to 20
  beta_parameters<-with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(het, list(m, c, gen), get_dist))
  beta_parameters<-cbind(beta_parameters[1:3], data.frame(beta_parameters[4]$x))
  names(beta_parameters)<-c("m", "c", "gen", "alpha", "beta")
  beta_parameters$index<-paste(beta_parameters$m, beta_parameters$c, beta_parameters$gen)
  alldata_df_6$index<-paste(alldata_df_6$m, alldata_df_6$c, alldata_df_6$gen)
  alldata_df_leftout<-alldata_df_6[alldata_df_6$rep==j, ]
  
  for(i in 1:length(beta_parameters$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dbeta(alldata_df_leftout[alldata_df_leftout$index==beta_parameters$index[i], ]$het, shape1=beta_parameters$alpha[i], shape2=beta_parameters$beta[i], log=TRUE))
  }
}

## give us a m, c, gen by replicate matrix. 

### I think that the next thing to do would be to make plots for everything, but all from one simulation at first? and then loop through all of the simulations

str(likelihoods)
summary(likelihoods)