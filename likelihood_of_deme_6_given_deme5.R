
### use the distribution of 19 deme 5s to predict 1 deme 6.
##### for number of junctions 
for(j in 1:20){ ##change this back to 20
  mean_njunct<-with(alldata_df[alldata_df$rep!=j & alldata_df$deme==5,], aggregate(njunct, list(m, c, gen, mech), mean))
  names(mean_njunct)<-c("m", "c", "gen", "mech", "mean_njunct")
  mean_njunct$index<-paste(mean_njunct$m, mean_njunct$c, mean_njunct$gen, mean_njunct$mech)#### there is definitely a cleverer way to do this
  alldata_df$index<-paste(alldata_df$m, alldata_df$c, alldata_df$gen, alldata_df$mech)
  alldata_df_leftout<-alldata_df[alldata_df$rep==j&alldata_df$deme==6, ]
  for(i in 1:length(mean_njunct$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dpois(alldata_df_leftout[alldata_df_leftout$index==mean_njunct$index[i], ]$njunct, lambda=mean_njunct$mean_njunct[i], log=TRUE))
     }
}

like_out <- cbind(mean_njunct[,1:4], likelihoods)
dim(like_out)
head(like_out)

###############
### for q #####
###############

### need to take the mean # of junctions for 1:19 replicates, for each generation, for each deme, for each m, for each c
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)

alldata_df$q<-ifelse(alldata_df$q ==0, 0.001, ifelse(alldata_df$q == 1, 0.999, alldata_df$q))


#### the current problem is that the starting shape really can't be the same for each of these. Which means that I can't run all 240 as a function, I've got to figure out what the starting shape is for each trial. 
get_dist<-function(x, shape1, shape2){
  estimates<-fitdistr(x, 'beta', start=list(shape1=shape1,shape2=shape2), lower=c(0.0001,0.0001)) ### not at all sure how to decide on the starting shapes!
  estimates$estimate
  ### this completely ignores the error around the alpha and beta distributions. Not sure about this at all. 
}

for(j in 1:20){ ##change this back to 20
  beta_parameters<-tryCatch(
    {
      message("This is the try part")
      with(alldata_df[alldata_df$rep!=j & alldata_df$deme==5,], aggregate(q, list(m, c, gen, mech), get_dist, shape1=5, shape2=5))
    }, 
    error=function(cond){
      message(paste("fitdistr isn't working because it doesn't like your starting values"))
      message("Here's the original error message:")
      message(cond)
      tryCatch(
        {
          message("This is the try again part")
          with(alldata_df[alldata_df$rep!=j & alldata_df$deme==5,], aggregate(q, list(m, c, gen, mech), get_dist, shape1=10, shape2=10))
        },
        error=function(again){
          message(paste("fitdistr isn't working because it still doesn't like your starting values"))
          message("Here's the original error message:")
          message(again)
          with(alldata_df[alldata_df$rep!=j & alldata_df$deme==5,], aggregate(q, list(m, c, gen, mech), get_dist, shape1=1, shape2=1))
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
  alldata_df$index<-paste(alldata_df$m, alldata_df$c, alldata_df$gen, alldata_df$mech)
  alldata_df_leftout<-alldata_df[alldata_df$rep==j&alldata_df$deme==6, ]
  
  for(i in 1:length(beta_parameters$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dbeta(alldata_df_leftout[alldata_df_leftout$index==beta_parameters$index[i], ]$q, shape1=beta_parameters$alpha[i], shape2=beta_parameters$beta[i], log=TRUE))
  }
}
