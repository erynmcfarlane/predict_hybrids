########################################
## 3) eryn's code for junctions
########################################

#10/02/22 attempt at asking how variable the number of junctions are across replicates of the same simulation.

#plan to extend this to q, Q and possibly individual fitness? 

require('stringr')

library(stringr)
library(MASS)
#library(fitdistrplus)
##if using local computer
datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/deme11_all" , pattern="*first8cols.txt.gz")
##if using teton
#datafiles<-list.files("/project/evolgen/jjahner/hybrid_sims/" , pattern="*first8cols.txt.gz")
m<-str_extract(datafiles, "(\\d+\\.*\\d*)")
c<-str_match(datafiles,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(datafiles, "^([^_]+_){1}([^_])") ### where e means there's an environmental interaction and m means there isn't

alldata<-list()

for(i in 1:length(datafiles)){
  alldata[[i]]<-read.table(gzfile(datafiles[i]), sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}
alldata_df<-do.call(rbind.data.frame, alldata)
summary(alldata_df)
### I'm not at all sure that I need this, so maybe take it out if it's not needed. 
alldata_df$gen<-as.factor(alldata_df$gen)
#need to take the mean # of junctions for 19 replicates, for each generation, for each deme, for each m, for each c, and ask how the individual number of junctions for the other replicate fits the distribution
#do this for each replicate, through all 20
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)
plots<-list()
#only going to use middle, hybridizing demes 5 and 6
#two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))
for(j in 1:20){ ##change this back to 20
  mean_njunct<-with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(njunct, list(m, c, gen, mech), mean))
  names(mean_njunct)<-c("m", "c", "gen", "deme", "mean_njunct")
  mean_njunct$index<-paste(mean_njunct$m, mean_njunct$c, mean_njunct$gen, mean_njunct$deme)#### there is definitely a cleverer way to do this
  alldata_df_6$index<-paste(alldata_df_6$m, alldata_df_6$c, alldata_df_6$gen, alldata_df_6$deme)
  alldata_df_leftout<-alldata_df_6[alldata_df_6$rep==j, ]
  for(i in 1:length(mean_njunct$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dpois(alldata_df_leftout[alldata_df_leftout$index==mean_njunct$index[i], ]$njunct, lambda=mean_njunct$mean_njunct[i], log=TRUE))
    ##I really don't want to do 600*20 plots.
    #plots[i]<-print(plot(dpois(alldata_df_leftout[alldata_df_leftout$index==mean_njunct$index[i], ]$njunct, lambda=mean_njunct$mean_njunct[i], log=TRUE)))
  }
}


like_out <- cbind(mean_njunct[,1:4], likelihoods)
dim(like_out)
head(like_out)

#write.table(like_out, file="deme11_junc_likes.txt", quote=F, row.names=F, col.names=F)




########################################
## 4) eryn's code for q
########################################


### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)

### only going to use middle, hybridizing demes 5 and 6
## two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))

### this is because the beta distribution is 0 to 1, excluding 0 and 1. I'm not entirely sure how much this matters, but I think it's what's needed to make fdistr beta to work
alldata_df_6$q<-ifelse(alldata_df_6$q ==0, 0.001, ifelse(alldata_df_6$q == 1, 0.999, alldata_df_6$q))


#### the current problem is that the starting shape really can't be the same for each of these. Which means that I can't run all 240 as a function, I've got to figure out what the starting shape is for each trial. 
get_dist<-function(x, shape1, shape2){

    estimates<-fitdistr(x, 'beta', start=list(shape1=shape1,shape2=shape2), lower=c(0.0001,0.0001)) ### not at all sure how to decide on the starting shapes!
  
  estimates$estimate
  ### this completely ignores the error around the alpha and beta distributions. Not sure about this at all. 
}

subset(alldata_df_6, mech=="dmi_m")->alldata_df_6
for(j in 1:20){ ##change this back to 20
  beta_parameters<-tryCatch(
    {
      message("This is the try part")
      with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(q, list(m, c, gen, mech), get_dist, shape1=5, shape2=5))
    }, 
    error=function(cond){
      message(paste("fitdistr isn't working because it doesn't like your starting values"))
      message("Here's the original error message:")
      message(cond)
      tryCatch(
        {
          message("This is the try again part")
          with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(q, list(m, c, gen, mech), get_dist, shape1=10, shape2=10))
        },
        error=function(again){
          message(paste("fitdistr isn't working because it still doesn't like your starting values"))
          message("Here's the original error message:")
          message(again)
          with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(q, list(m, c, gen, mech), get_dist, shape1=1, shape2=1))
        }, 
        warning=function(again){
          message(paste("You've got a warning, that's fine"))
          message("Here's the warning:")
          message(again)
        }, 
        finally={
          message(paste("It's over"))
        }
      )
    },
    warning=function(cond){
      message(paste("You've got a warning, that's fine"))
      message("Here's the warning:")
      message(cond)
    }, 
    finally={
      message(paste("It's over"))
    })
  beta_parameters<-cbind(beta_parameters[1:4], data.frame(beta_parameters[5]$x))
  names(beta_parameters)<-c("m", "c", "gen","mechanism","alpha", "beta")
  beta_parameters$index<-paste(beta_parameters$m, beta_parameters$c, beta_parameters$gen)
  alldata_df_6$index<-paste(alldata_df_6$m, alldata_df_6$c, alldata_df_6$gen, alldata_df_6$mech)
  alldata_df_leftout<-alldata_df_6[alldata_df_6$rep==j, ]
  
  for(i in 1:length(beta_parameters$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dbeta(alldata_df_leftout[alldata_df_leftout$index==beta_parameters$index[i], ]$q, shape1=beta_parameters$alpha[i], shape2=beta_parameters$beta[i], log=TRUE))
  }
}

## give us a m, c, gen by replicate matrix. 

### I think that the next thing to do would be to make plots for everything, but all from one simulation at first? and then loop through all of the simulations

str(likelihoods)
summary(likelihoods)

likes <- cbind(beta_parameters[1:4], likelihoods)
dim(likes)
head(likes)

#write.table(likes, file="deme11_q_likes.txt", quote=F, row.names=F, col.names=F)

########################################
## 5) eryn's code for het (Q)
########################################

### doing the same thing, but for het, using a beta. Many, many more questions. 

### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)

### only going to use middle, hybridizing demes 5 and 6
## two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))

### this is because the beta distribution is 0 to 1, excluding 0 and 1. I'm not entirely sure how much this matters, but I think it's what's needed to make fdistr beta to work
alldata_df_6$het<-ifelse(alldata_df_6$het ==0, alldata_df_6$het+0.001, ifelse(alldata_df_6$het == 1, alldata_df_6$het-0.001, alldata_df_6$het))

for(j in 1:20){ ##change this back to 20
  beta_parameters<-tryCatch(
    {
      message("This is the try part")
      with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(het, list(m, c, gen, mech), get_dist, shape1=5, shape2=5))
    }, 
    error=function(cond){
      message(paste("fitdistr isn't working because it doesn't like your starting values"))
      message("Here's the original error message:")
      message(cond)
      tryCatch(
        {
          message("This is the try again part")
          with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(het, list(m, c, gen, mech), get_dist, shape1=10, shape2=10))
        },
        error=function(again){
          message(paste("fitdistr isn't working because it still doesn't like your starting values"))
          message("Here's the original error message:")
          message(again)
          with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(het, list(m, c, gen, mech), get_dist, shape1=1, shape2=1))
        }, 
        warning=function(again){
          message(paste("You've got a warning, that's fine"))
          message("Here's the warning:")
          message(again)
        }, 
        finally={
          message(paste("It's over"))
        }
      )
    },
    warning=function(cond){
      message(paste("You've got a warning, that's fine"))
      message("Here's the warning:")
      message(cond)
    }, 
    finally={
      message(paste("It's over"))
    })
  beta_parameters<-cbind(beta_parameters[1:4], data.frame(beta_parameters[5]$x))
  names(beta_parameters)<-c("m", "c", "gen","mechanism","alpha", "beta")
  beta_parameters$index<-paste(beta_parameters$m, beta_parameters$c, beta_parameters$gen)
  alldata_df_6$index<-paste(alldata_df_6$m, alldata_df_6$c, alldata_df_6$gen, alldata_df_6$mech)
  alldata_df_leftout<-alldata_df_6[alldata_df_6$rep==j, ]
  
  for(i in 1:length(beta_parameters$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dbeta(alldata_df_leftout[alldata_df_leftout$index==beta_parameters$index[i], ]$het, shape1=beta_parameters$alpha[i], shape2=beta_parameters$beta[i], log=TRUE))
  }
}

## give us a m, c, gen by replicate matrix. 

### I think that the next thing to do would be to make plots for everything, but all from one simulation at first? and then loop through all of the simulations

str(likelihoods)
summary(likelihoods)

likes <- cbind(beta_parameters[1:3], likelihoods)
dim(likes)
head(likes)

#write.table(likes, file="deme11_het_likes.txt", quote=F, row.names=F, col.names=F)



########################################
## 6) plotting likelihoods
########################################

#likes <- read.delim("deme11_junc_likes.txt", sep=" ", header=FALSE)
#likes <- read.delim("deme11_q_likes.txt", sep=" ", header=FALSE)
likes <- read.delim("deme11_het_likes.txt", sep=" ", header=FALSE)
dim(likes)
head(likes)


uniq_m <- unique(likes[,1])
uniq_c <- unique(likes[,2])
uniq_gen <- unique(likes[,3])

quartz(height=7.5, width=11.25)
par(mar=c(5,5,2,2), mfrow=c(2,3))


for (j in 1:length(uniq_m))
{
  m_sub <- subset(likes, likes[,1]==uniq_m[j])
  for (k in 1:length(uniq_c))
  {
    c_sub <- subset(m_sub, m_sub[,2]==uniq_c[k])
    plot(0, type="n", xlim=c(0,length(uniq_gen)), ylim=c(range(likes[,4:23])[1],range(likes[,4:23])[2]), xlab="Generation", ylab="Joint log likelihood", cex.lab=1.5, xaxt="n")
    axis(1, at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5), labels=uniq_gen)
    mtext(paste0("m: ", uniq_m[j], "; c: ", uniq_c[k]), line=0.5)
    box(lwd=2)
    for (l in 1:length(uniq_gen))
    {
      gen_sub <- subset(c_sub, c_sub[,3]==uniq_gen[l])
      #segments(l-0.5, mean(as.numeric(gen_sub[1,5:24]))-sd(as.numeric(gen_sub[1,5:24])), l-0.5, mean(as.numeric(gen_sub[1,5:24]))+sd(as.numeric(gen_sub[1,5:24])), lwd=1.5)
      #points(l-0.5, mean(as.numeric(gen_sub[1,5:24])), pch=21, bg="gray", cex=2, lwd=1.25)
      for (m in 4:23) { points(l-0.5, gen_sub[1,m], pch=21, bg=adjustcolor("gray", alpha.f=0.5), cex=1.5) }
    }
  }
}
