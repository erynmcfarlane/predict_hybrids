#### interpretation of the genotypes likelihood script####
load("Predicting_hybrid_genotypes.RData")
### sense check - did it work? ###
summary(likelihoods) #### therea re sooooo many NAs

rownames(likelihoods)<-mean_freq$index

likelihoods_df[likelihoods_df[,1]==max(likelihoods_df[,1]),]
likelihoods_df[likelihoods_df[,2]==max(likelihoods_df[,2]),]
likelihoods_df[likelihoods_df[,3]==max(likelihoods_df[,3]),]

### which SNPs have the DMI?###