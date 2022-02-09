### first attempt at random forests for predicting hybrid scores

### question: can we predict individual Q and q scores (which is what we can measure empirically) using model parameters from Lindkt and Buerkle 2015?
### essentially want to ask if the Q and q scores are reflective of the processes that we know are happening, given that the simulations stipulate migration, selection against hybrids and genetic mechanism leading to hybrid dysfunction

###when working on local computer
setwd("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/predict_hybrids/")

##########################
#### load in the data#####
##########################

datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/" , pattern="*first8cols.txt.gz")

alldata<-list()
for(i in 1:length(datafiles)){
  alldata[[i]]<-read.table(gzfile(datafiles[i]), sep=",", header=T)
}

alldata_df<-do.call(rbind.data.frame, alldata)

#######################################
####load packages for random forest####
#######################################


