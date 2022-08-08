### first attempt at random forests for predicting hybrid scores

#authors include Eryn McFarlane, Josh Jahner, Alex Buerkle, Liz Mandeville

### question: can we predict individual Q and q scores (which is what we can measure empirically) using model parameters from Lindkt and Buerkle 2015?
### essentially want to ask if the Q and q scores are reflective of the processes that we know are happening, given that the simulations stipulate migration, selection against hybrids and genetic mechanism leading to hybrid dysfunction

###when working on local computer
setwd("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/")

##########################
#### load in the data#####
##########################
library(stringr)
datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/" , pattern="*first8cols.txt.gz")
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
### at the moment, I'm going to subset the data, because I can't run this. 

subset_data<-alldata_df[sample(nrow(alldata_df), nrow(alldata_df)*0.05),]
#######################################
####load packages for random forest####
#######################################

library(randomForest)
library(varSelRF)
library(randomForestExplainer)
library(ROCR)
library(caret)

### we can use 90% (or something) of the data to predict the other 10%
### alternatively, we can do a 10-fold cross validation, and do this 10 times (9 times?)

#### use 90% of the data to predict the other 10% of the data####
nrow(subset_data)

training_pop<-subset_data[sample(nrow(subset_data), nrow(subset_data)*0.9),]
testing_pop<-subset(subset_data, ! (rownames(subset_data) %in% rownames(training_pop)))

rf_hybrids <- randomForest(
  x = training_pop[,c(1:4, 7:10)], 							#these are the predictors ## leave out variable 1
  y = training_pop[,5]),						#response - I've made it factor here to force classification
  ### I made the c a numeric above because I've not checked how variable it is over all
  importance = T,						#save variable importance metrics
  ntree=400,							#number of trees to make - general rules, go for 10x's the number of features
  nodesize = 5, ###performance can change with different node sizes
  mtry=5)	    

rf_hybrids

### want to maximize the Gini Gain. 
varImpPlot(rf_hybrids, sort=T)


out <- predict(rf_hybrids,testing_pop,'response') #this is the syntax you would use if you wanted to predict to new data..
#just put new data frame in place of "test"..response just stays response so that it knows to generate yhat

plot(out~testing_pop$q, xlab="real q", ylab="predicted q")

#### I'm a little worried that I have the plan backwards. So I'm going to see if I can use q,Q, junctions and deme to predict something.

rf_hybrids_predict_c <- randomForest(
  x = training_pop[,1:9], 							#these are the predictors
  y = as.factor(training_pop[,10]),						##response - I've made it factor here to force classification
  ### I made the c a numeric above because I've not checked how variable it is over all
  importance = T,						#save variable importance metrics
  ntree=400,							#number of trees to make - general rules, go for 10x's the number of features
  nodesize = 9, ###performance can change with different node sizes ### this definitely needs to be played around with, given how much data there are. 
  mtry=5)	    

rf_hybrids_predict_c ### this gives a confusion matrix

#### the error rate is really high!!!

varImpPlot(rf_hybrids_predict_c, sort=T)

out <- predict(rf_hybrids_predict_c,testing_pop,'response')
### I'm not entirely sure what this is telling me, or how to interpret the confusion matrix
plot(out~testing_pop$c, xlab="real c", ylab="predicted c")



####################################################################################
##### now that I've figured out the jist of how to do this, I want to tune these####
####################################################################################
### following this https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/
### start tuning mtry - number of variables randomly sampled as candidates
### and ntree - number of trees to grow - obs this makes a huge difference on how long this takes to run
set.seed(42)
y<-training_pop[,10]
x<-training_pop[,1:9]
bestmtry <- tuneRF(x, y, stepFactor=1.5, improve=1e-5, ntree=500)
print(bestmtry)

