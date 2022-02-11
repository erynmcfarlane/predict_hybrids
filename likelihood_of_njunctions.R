dpois(c(1,2,5,12), lambda=5, log=TRUE)
#this gives the point probabilities of each othese in the poisson distriubtions
sum(dpois(c(1,2,5,12), lambda=5, log=TRUE))

### take the mean of the variable of interest for the lambda, and then take all of the individuals for the 20th replicate to then see
### for this replicate, what's the likelihood of these individuals
### if there's high variablity, the liklihood is gonna suck.
### this is gonna vary depending on the simulations
## we could use this as a measure of variablity. 
### we can do this for junctions, intersouce ancestry and little q

### intersource ancesty and litle q will need a beta distribution,
### beta(alpha, beta), mean=alpha/(alpha+ beta)
### beta(pi, theta)

#beta(pi * theta, (1-pi) * theta))
#mean=pi
#precision is theta, theta = 1/var

### can use fitdist() with beta to characterize the 19 replicates

### time for a loop, get a distribution of joint likihoods. 
### get a distribution of log likelihoods where we've left it out
### contrast this accross time, migration, selection
### only do this for one deme, the middle deme, for deme 11, where there is a middle deme (or deme 3)



#### the plan is to see if we can use 19 demes from a simulation to predict the 20th, accross different generations. 
###when working on local computer
setwd("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/")

##########################
#### load in the data#####
##########################
library(stringr)
##if using local computer
datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/" , pattern="*first8cols.txt.gz")
##if using teton
#datafiles<-list.files("/project/evolgen/jjahner/hybrid_sims/" , pattern="*first8cols.txt.gz")
m<-str_extract(datafiles, "(\\d+\\.*\\d*)")
c<-str_match(datafiles,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?

alldata<-list()
for(i in 1:length(datafiles)){
  alldata[[i]]<-read.table(gzfile(datafiles[i]), sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]])))
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
}

alldata_df<-do.call(rbind.data.frame, alldata)

summary(alldata_df)

### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c

likelihoods<-matrix(nrow=120, ncol=20)
plots<-list()

### only going to use middle, hybridizing demes 5 and 6
## two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_56<-subset(alldata_df, deme %in% c(5,6))

for(j in 1:20){ ##change this back to 20
mean_njunct<-with(alldata_df_56[alldata_df_56$rep!=j,], aggregate(njunct, list(m, c, gen, deme), mean))
names(mean_njunct)<-c("m", "c", "gen", "deme", "mean_njunct")
mean_njunct$index<-paste(mean_njunct$m, mean_njunct$c, mean_njunct$gen, mean_njunct$deme)#### there is definitely a cleverer way to do this
alldata_df_56$index<-paste(alldata_df_56$m, alldata_df_56$c, alldata_df_56$gen, alldata_df_56$deme)
alldata_df_leftout<-alldata_df_56[alldata_df_56$rep==j, ]

for(i in 1:length(mean_njunct$index)){ ##change this back to length(mean_njunct$index)
likelihoods[i,j]<-sum(dpois(alldata_df_leftout[alldata_df_leftout$index==mean_njunct$index[i], ]$njunct, lambda=mean_njunct$mean_njunct[i], log=TRUE))
##I really don't want to do 600*20 plots.
#plots[i]<-print(plot(dpois(alldata_df_leftout[alldata_df_leftout$index==mean_njunct$index[i], ]$njunct, lambda=mean_njunct$mean_njunct[i], log=TRUE)))
}
}

## give us a m, c, gen, deme by replicate matrix. 

### I think that the next thing to do would be to make plots for everything, but all from one simulation at first? and then loop through all of the simulations

str(likelihoods)
