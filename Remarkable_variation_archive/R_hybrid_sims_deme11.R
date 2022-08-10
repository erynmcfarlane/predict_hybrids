########################################
## R_hybrid_sims.R
########################################

## JPJ 2 ii 22

## What's in this file?
	## 1) load data and packages
	## 2) initial plotting exercise
	## 3) eryn's code for junctions
	## 4) eryn's code for q
	## 5) eryn's code for Q
	## 6) plotting likelihoods for junctions (old)
	## 7) plotting distribution case studies
	## 8) predicting future generations
	## 9) density figures
	## 10) pca
	## 11) anova plots
	## 12) contingency
	## 13) qQ across replicates
	## 14) gen10 vs gen100


########################################
## 1) load data and packages
########################################

library(MetBrewer)


########################################
## 2) initial plotting exercise
########################################

sim_data <- read.csv("dmi_m0.01_c0.9.main_first8cols.txt", header=TRUE)
	dim(sim_data)
	head(sim_data)

rep_list <- unique(sim_data[,1])
gen_list <- unique(sim_data[,2])
deme_list <- rev(unique(sim_data[,3]))
ind_list <- unique(sim_data[,4])


for (i in 1:length(rep_list))
	{
	pdf(paste0("dmi_m0.2_neutral_qQ_rep",i,".pdf"), width=25, height=25)
	par(mar=c(5,5,1,1), mfrow=c(10,10))
	rep_sub <- subset(sim_data, sim_data[,1]==rep_list[i])
	for (j in 1:length(gen_list))
		{
		gen_sub <- subset(rep_sub, rep_sub[,2]==gen_list[j])
		for (k in 1:length(deme_list))
			{
			deme_sub <- subset(gen_sub, gen_sub[,3]==deme_list[k])
			plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="q", ylab="het", las=1, cex.lab=1.5, cex.axis=1.25)
			segments(0, 0, 0.5, 1, lty=2, lwd=3)
			segments(1, 0, 0.5, 1, lty=2, lwd=3)
			mtext(paste0("Rep: ", i), cex=0.75, adj=0.02, line=-1.5)
			mtext(paste0("Gen: ", j), cex=0.75, adj=0.02, line=-3)
			mtext(paste0("Deme: ", k), cex=0.75, adj=0.02, line=-4.5)
			box(lwd=2)
			for (l in 1:length(ind_list))
				{
				points(deme_sub[l,5], deme_sub[l,6], pch=21, bg=adjustcolor("red", alpha.f=0.5), cex=2)
				}
			print(c(i,j,k))
			}
		}
	dev.off()
	}



########################################
## 3) eryn's code for junctions
########################################

#10/02/22 attempt at asking how variable the number of junctions are across replicates of the same simulation.

#plan to extend this to q, Q and possibly individual fitness? 

require('stringr')

#load in the data from teton


datafiles<-list.files("/project/evolgen/jjahner/hybrid_sims/deme11_3/" , pattern="first8cols.txt.gz")
m<-str_extract(datafiles, "(\\d+\\.\\d*)")
c<-str_match(datafiles,"c(\\d+\\.\\d)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?

alldata<-list()

for(i in 1:length(datafiles)){
alldata[[i]]<-read.table(gzfile(datafiles[i]), sep=",", header=T)
alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
}
alldata_df<-do.call(rbind.data.frame, alldata)
summary(alldata_df)
#need to take the mean # of junctions for 19 replicates, for each generation, for each deme, for each m, for each c, and ask how the individual number of junctions for the other replicate fits the distribution
#do this for each replicate, through all 20
likelihoods<-matrix(nrow=60, ncol=20)
plots<-list()
#only going to use middle, hybridizing demes 5 and 6
#two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))
for(j in 1:20){ ##change this back to 20
mean_njunct<-with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(njunct, list(m, c, gen, deme), mean))
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

like_out <- cbind(mean_njunct[,1:3], likelihoods)
	dim(like_out)
	head(like_out)
	
#write.table(like_out, file="deme11_3_junc_likes.txt", quote=F, row.names=F, col.names=F)




########################################
## 4) eryn's code for q
########################################

### doing the same thing, but for q, using a beta. Many, many more questions. 

##########################
#### load in the data#####
##########################
library(stringr)
library(MASS)
library(fitdistrplus)
##if using local computer
#datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/deme11" , pattern="*first8cols.txt.gz")
##if using teton - something needs to change here to get data from all 3 folders. 
datafiles<-list.files("/project/evolgen/jjahner/hybrid_sims/deme11_3/" , pattern="*first8cols.txt.gz")
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
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)

### only going to use middle, hybridizing demes 5 and 6
## two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))

### this is because the beta distribution is 0 to 1, excluding 0 and 1. I'm not entirely sure how much this matters, but I think it's what's needed to make fdistr beta to work
alldata_df_6$q<-ifelse(alldata_df_6$q ==0, alldata_df_6$q+0.001, ifelse(alldata_df_6$q == 1, alldata_df_6$q-0.001, alldata_df_6$q))

get_dist<-function(x){
  estimates<-fitdistr(x, 'beta', start=list(shape1=10,shape2=10)) ### not at all sure how to decide on the starting shapes!
  estimates$estimate
  ### this completely ignores the error around the alpha and beta distributions. Not sure about this at all. 
}

for(j in 1:20){ ##change this back to 20
  beta_parameters<-with(alldata_df_6[alldata_df_6$rep!=j,], aggregate(q, list(m, c, gen), get_dist))
  beta_parameters<-cbind(beta_parameters[1:3], data.frame(beta_parameters[4]$x))
  names(beta_parameters)<-c("m", "c", "gen", "alpha", "beta")
  beta_parameters$index<-paste(beta_parameters$m, beta_parameters$c, beta_parameters$gen)
  alldata_df_6$index<-paste(alldata_df_6$m, alldata_df_6$c, alldata_df_6$gen)
  alldata_df_leftout<-alldata_df_6[alldata_df_6$rep==j, ]
  
  for(i in 1:length(beta_parameters$index)){ ##change this back to length(mean_njunct$index)
    likelihoods[i,j]<-sum(dbeta(alldata_df_leftout[alldata_df_leftout$index==beta_parameters$index[i], ]$q, shape1=beta_parameters$alpha[i], shape2=beta_parameters$beta[i], log=TRUE))
  }
}

## give us a m, c, gen by replicate matrix. 

### I think that the next thing to do would be to make plots for everything, but all from one simulation at first? and then loop through all of the simulations

str(likelihoods)
summary(likelihoods)

likes <- cbind(beta_parameters[1:3], likelihoods)
	dim(likes)
	head(likes)

#write.table(likes, file="deme11_3_q_likes.txt", quote=F, row.names=F, col.names=F)




########################################
## 5) eryn's code for het (Q)
########################################

### doing the same thing, but for het, using a beta. Many, many more questions. 

##########################
#### load in the data#####
##########################
library(stringr)
library(MASS)
library(fitdistrplus)
##if using local computer
#datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/deme11" , pattern="*first8cols.txt.gz")
##if using teton - something needs to change here to get data from all 3 folders. 
datafiles<-list.files("/project/evolgen/jjahner/hybrid_sims/deme11_3/" , pattern="*first8cols.txt.gz")
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
nrow<-10*length(datafiles)

likelihoods<-matrix(nrow=nrow, ncol=20)

### only going to use middle, hybridizing demes 5 and 6
## two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))

### this is because the beta distribution is 0 to 1, excluding 0 and 1. I'm not entirely sure how much this matters, but I think it's what's needed to make fdistr beta to work
alldata_df_6$het<-ifelse(alldata_df_6$het ==0, alldata_df_6$het+0.001, ifelse(alldata_df_6$het == 1, alldata_df_6$het-0.001, alldata_df_6$het))

get_dist<-function(x){
  estimates<-fitdistr(x, 'beta', start=list(shape1=10,shape2=10)) ### not at all sure how to decide on the starting shapes!
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

likes <- cbind(beta_parameters[1:3], likelihoods)
	dim(likes)
	head(likes)

#write.table(likes, file="deme11_3_het_likes.txt", quote=F, row.names=F, col.names=F)




########################################
## 6) plotting likelihoods
########################################

likes <- read.delim("deme11_3_junc_likes.txt", sep=" ", header=FALSE)
#likes <- read.delim("deme11_3_q_likes.txt", sep=" ", header=FALSE)
#likes <- read.delim("deme11_3_het_likes.txt", sep=" ", header=FALSE)
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




########################################
## 7) plotting histograms
########################################
	
plot_hists(trial_in="dmi", m_in=0.01, c_in=0.9, gen_in=100, rep_in=6)
	## input file for example: dmi_m0.01_c0.9.main_first8cols.txt
	
plot_hists <- function(trial_in=NA, m_in=NA, c_in=NA, gen_in=NA, rep_in=NA)
	## potential trial_in values: dmi, path, dmi_e, path_e
	## potential m_in values: 0.01, 0.2
	## potential c_in values: neutral, 0.2, 0.9
	## potential gen_in values: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100
	## potential rep_in values: 1:20
	{
	require(fitdistrplus)
	require(plyr)
	require(MetBrewer)
	hist_cols <- met.brewer("Derain", 7)

	## create file handle
	if (c_in=="neutral")		{ handle <- paste0(trial_in, "_m", m_in, "_", c_in, ".main_first8cols.txt") }
	else						{ handle <- paste0(trial_in, "_m", m_in, "_c", c_in, ".main_first8cols.txt") }
	
	## load and subset data
	dat_in <- read.csv(handle, header=TRUE)
	deme_sub <- subset(dat_in, dat_in[,3]==6)
	gen_sub <- subset(deme_sub, deme_sub[,2]==gen_in)
	gen_sub_rep <- subset(gen_sub, gen_sub[,1]==rep_in)
	gen_sub_others <- subset(gen_sub, gen_sub[,1]!=rep_in)
	
	## generate beta_distribution lines based on the 19 other replicates
	q_beta <-fitdistr(gen_sub_others[,5], 'beta', start=list(shape1=10,shape2=10))
	Q_beta <-fitdistr(gen_sub_others[,6], 'beta', start=list(shape1=10,shape2=10))
	#junc_negbin <-fitdistr(gen_sub_others[,8], 'negative binomial')

	## start plot
	if (c_in=="neutral")		{ pdf(paste0("hist_", trial_in, "_m", m_in, "_", c_in, "_gen", gen_in, "_rep", rep_in, ".pdf"), height=10, width=10) }
	else						{ pdf(paste0("hist_", trial_in, "_m", m_in, "_c", c_in, "_gen", gen_in, "_rep", rep_in, ".pdf"), height=10, width=10) }
	#quartz(height=10, width=10)
	par(mar=c(5,5,1,1), mfrow=c(2,2))
	
	## junction histogram
	upper_break <- round_any(max(gen_sub[,8]), 10, ceiling)
	lower_break <- round_any(min(gen_sub[,8]), 10, floor)
	break_vect_junc <- seq(lower_break, upper_break, 10)
	plot(seq(0, upper_break, 1), dpois(seq(0, upper_break, 1), mean(gen_sub_others[,8])), type="l", axes=FALSE, xlab="", ylab="", lwd=2, lty=2)
	par(new=TRUE)
	plot(0, type="n", xlab="Number of junctions", ylab="Number of individuals", xlim=c(lower_break, upper_break), ylim=c(0, 150), las=1, cex.lab=1.5, cex.axis=1.25)
	hist(gen_sub_others[,8], breaks=break_vect_junc, col=adjustcolor(hist_cols[5], alpha.f=0.7), add=TRUE, lwd=1.25)
	hist(gen_sub_rep[,8], breaks= break_vect_junc, col=adjustcolor(hist_cols[1], alpha.f=0.7), add=TRUE, lwd=1.25)
	mtext("A", cex=2, adj=-0.2, line=-1)
	box(lwd=2)


	## qQ plot
	plot(0, type="n", xlab=expression("Admixture proportion ("~italic("q")~")"), ylab=expression("Inter-source ancestry ("~italic("Q")[12]~")"), xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
	segments(0, 0, 0.5, 1, lwd=4, lty=2)
	segments(1, 0, 0.5, 1, lwd=4, lty=2)
	points(gen_sub_others[,5], gen_sub_others[,6], pch=21, bg=adjustcolor(hist_cols[5], alpha.f=0.7), cex=1.5)
	points(gen_sub_rep[,5], gen_sub_rep[,6], pch=21, bg=adjustcolor(hist_cols[1], alpha.f=0.7), cex=1.5)
	legend("topleft", legend=c(paste0("Replicate ", rep_in), "Other replicates"), pch=22, pt.bg=c(adjustcolor(hist_cols[1], alpha.f=0.7), adjustcolor(hist_cols[5], alpha.f=0.7)), cex=1.25, pt.cex=2.5, bty="n")
	mtext(paste0("trial: ", trial_in, "; m: ", m_in), line=-1.9, adj=0.98)
	mtext(paste0("c: ", c_in, "; gen: ", gen_in), line=-3.15, adj=0.98)
	mtext("B", cex=2, adj=-0.2, line=-1)
	box(lwd=2)
	
	## q histogram
	plot(seq(0,1,0.01), dbeta(seq(0,1,0.01), q_beta$estimate[1], q_beta$estimate[2]), type="l", axes=FALSE, xlab="", ylab="", lwd=2, lty=2)
	par(new=TRUE)
	break_vect <- seq(0, 1, by=0.01)
	plot(0, type="n", xlab=expression("Admixture proportion ("~italic("q")~")"), ylab="Number of individuals", xlim=c(0, 1), ylim=c(0, 200), las=1, cex.lab=1.5, cex.axis=1.25)
	hist(gen_sub_others[,5], breaks=break_vect, col=adjustcolor(hist_cols[5], alpha.f=0.7), add=TRUE, lwd=1.25)
	hist(gen_sub_rep[,5], breaks=break_vect, col=adjustcolor(hist_cols[1], alpha.f=0.7), add=TRUE, lwd=1.25)
	mtext("C", cex=2, adj=-0.2, line=-1)
	box(lwd=2)
	
	## het histogram
	plot(seq(0,1,0.01), dbeta(seq(0,1,0.01), Q_beta$estimate[1], Q_beta$estimate[2]), type="l", axes=FALSE, xlab="", ylab="", lwd=2, lty=2)
	par(new=TRUE)
	plot(0, type="n", xlab=expression("Inter-source ancestry ("~italic("Q")[12]~")"), ylab="Number of individuals", xlim=c(0, 1), ylim=c(0, 200), las=1, cex.lab=1.5, cex.axis=1.25)
	hist(gen_sub_others[,6], breaks=break_vect, col=adjustcolor(hist_cols[5], alpha.f=0.7), add=TRUE, lwd=1.25)
	hist(gen_sub_rep[,6], breaks=break_vect, col=adjustcolor(hist_cols[1], alpha.f=0.7), add=TRUE, lwd=1.25)
	mtext("D", cex=2, adj=-0.2, line=-1)
	box(lwd=2)
	dev.off()
	
	## summary stats to print
	rep_mean_junc <- round(mean(gen_sub_rep[,8]), digits=2)
	other_mean_junc <- round(mean(gen_sub_others[,8]), digits=2)
	rep_mean_q <- round(mean(gen_sub_rep[,5]), digits=2)
	other_mean_q <- round(mean(gen_sub_others[,5]), digits=2)
	rep_mean_het <- round(mean(gen_sub_rep[,6]), digits=2)
	other_mean_het <- round(mean(gen_sub_others[,6]), digits=2)
	print(paste0("trial: ", trial_in, "; m: ", m_in, "; c: ", c_in, "; gen: ", gen_in))
	print(paste0("rep ", rep_in, " mean junctions = ", rep_mean_junc, "; other reps mean junctions = ", other_mean_junc))
	print(paste0("rep ", rep_in, " mean q = ", rep_mean_q, "; other reps mean q = ", other_mean_q))
	print(paste0("rep ", rep_in, " mean het = ", rep_mean_het, "; other reps mean het = ", other_mean_het))
	}




########################################
## 8) predicting future generations
########################################

require(fitdistrplus)
require(plyr)
require(MetBrewer)


uniq_trials <- c("dmi", "path", "dmi_e", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c("neutral", 0.2, 0.9)
uniq_gen <- seq(20, 100, 10)
uniq_rep <- 1:20
uniq_metric <- c("njunct", "q", "het")


#for (i in 1:length(uniq_trials))
for (i in 1:1)
	{
	quartz(height=18, width=9)
	par(mar=c(5,5,1,1), mfrow=c(6,3))
	for (j in 1:length(uniq_m))
		{
		for (k in 1:length(uniq_c))
			{
			if (uniq_c[k]=="neutral")	{ handle <- paste0(uniq_trials[i], "_m", uniq_m[j], "_", uniq_c[k], ".main_first8cols.txt") }
			else							{ handle <- paste0(uniq_trials[i], "_m", uniq_m[j], "_c", uniq_c[k], ".main_first8cols.txt") }
			dat_in <- read.csv(handle, header=TRUE)
			deme_sub <- subset(dat_in, dat_in[,3]==6)
			for (z in 1:dim(deme_sub))
				{
				if (deme_sub[z,5]==0) { deme_sub[z,5] <- 0.001 }
				else if (deme_sub[z,5]==1) { deme_sub[z,5] <- 0.999 }
				}
			for (z in 1:dim(deme_sub))
				{
				if (deme_sub[z,6]==0) { deme_sub[z,6] <- 0.001 }
				else if (deme_sub[z,6]==1) { deme_sub[z,6] <- 0.999 }
				}
			gen_10 <- subset(deme_sub, deme_sub[,2]==10)
			gen_others <- subset(deme_sub, deme_sub[,2]!=10)
			for (l in 1:length(uniq_metric))
				{
				if (l==1)
					{
					plot(0, type="n", xlim=c(20,100), ylim=c(-200000, 0), xlab="Generation", ylab="Junctions JLL", cex.lab=1.5, xaxt="n")
					axis(1, at=seq(20,100,10), labels=seq(20,100,10))
					for (m in 1:length(uniq_rep))
						{
						gen_10_rep_sub <- subset(gen_10, gen_10[,1]==uniq_rep[m])
						gen_others_rep_sub <- subset(gen_others, gen_others[,1]==uniq_rep[m])
						for (n in 1:length(uniq_gen))
							{
							gen_others_rep_sub_gen_sub <- subset(gen_others_rep_sub, gen_others_rep_sub[,2]==uniq_gen[n])
							junc_lik <- sum(dpois(gen_others_rep_sub_gen_sub[,8], lambda=mean(gen_10_rep_sub[,8]), log=TRUE))
							points(uniq_gen[n], junc_lik, pch=21, bg=adjustcolor("grey", alpha.f=0.75), cex=2)
							}
						}
					}
				else if (l==2)
					{
					plot(0, type="n", xlim=c(20,100), ylim=c(-200, 500), xlab="Generation", ylab="q JLL", cex.lab=1.5, xaxt="n")
					axis(1, at=seq(20,100,10), labels=seq(20,100,10))
					for (m in 1:length(uniq_rep))
						{
						gen_10_rep_sub <- subset(gen_10, gen_10[,1]==uniq_rep[m])
						gen_others_rep_sub <- subset(gen_others, gen_others[,1]==uniq_rep[m])
						for (n in 1:length(uniq_gen))
							{
							gen_others_rep_sub_gen_sub <- subset(gen_others_rep_sub, gen_others_rep_sub[,2]==uniq_gen[n])
							q_fit <- fitdistr(gen_10_rep_sub[,5], 'beta', start=list(shape1=10,shape2=10))
							q_lik <- sum(dbeta(gen_others_rep_sub_gen_sub[,5], shape1=q_fit$estimate[1], shape2=q_fit$estimate[2], log=TRUE))
							points(uniq_gen[n], q_lik, pch=21, bg=adjustcolor("grey", alpha.f=0.75), cex=2)
							}
						}
					}
				else if (l==3)
					{
					plot(0, type="n", xlim=c(20,100), ylim=c(-200, 500), xlab="Generation", ylab="Q JLL", cex.lab=1.5, xaxt="n")
					axis(1, at=seq(20,100,10), labels=seq(20,100,10))
					mtext(paste0("m: ", uniq_m[j], "; c: ", uniq_c[k]), line=-1.9, adj=0.98)
					for (m in 1:length(uniq_rep))
						{
						gen_10_rep_sub <- subset(gen_10, gen_10[,1]==uniq_rep[m])
						gen_others_rep_sub <- subset(gen_others, gen_others[,1]==uniq_rep[m])
						for (n in 1:length(uniq_gen))
							{
							gen_others_rep_sub_gen_sub <- subset(gen_others_rep_sub, gen_others_rep_sub[,2]==uniq_gen[n])
							Q_fit <- fitdistr(gen_10_rep_sub[,6], 'beta', start=list(shape1=10,shape2=10))
							Q_lik <- sum(dbeta(gen_others_rep_sub_gen_sub[,6], shape1=Q_fit$estimate[1], shape2=Q_fit$estimate[2], log=TRUE))
							points(uniq_gen[n], Q_lik, pch=21, bg=adjustcolor("grey", alpha.f=0.75), cex=2)
							}
						}
					}
				}
			}
		}
	}
	
	



########################################
## 9) density figures
########################################

	
plot_densities(trial_in="dmi", m_in=0.01, c_in=0.9, gen_in=100)
	## input file for example: dmi_m0.01_c0.9.main_first8cols.txt
	
plot_densities <- function(trial_in=NA, m_in=NA, c_in=NA, gen_in=NA)
	## potential trial_in values: dmi, path, dmi_e, path_e
	## potential m_in values: 0.01, 0.2
	## potential c_in values: neutral, 0.2, 0.9
	## potential gen_in values: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100
	{

	## create file handle
	if (c_in=="neutral")		{ handle <- paste0(trial_in, "_m", m_in, "_", c_in, ".main_first8cols.txt") }
	else						{ handle <- paste0(trial_in, "_m", m_in, "_c", c_in, ".main_first8cols.txt") }
	
	## load and subset data
	dat_in <- read.csv(handle, header=TRUE)
	deme_sub <- subset(dat_in, dat_in[,3]==6)
	gen_sub <- subset(deme_sub, deme_sub[,2]==gen_in)
	
	## start plot
	if (c_in=="neutral")		{ pdf(paste0("density_", trial_in, "_m", m_in, "_", c_in, "_gen", gen_in, ".pdf"), height=4, width=12) }
	else						{ pdf(paste0("density_", trial_in, "_m", m_in, "_c", c_in, "_gen", gen_in, ".pdf"), height=4, width=12) }
	#quartz(height=4, width=12)
	par(mar=c(5,5,1,1), mfrow=c(1,3))
	
	## junction densities
	plot(0, type="n", xlab="Number of junctions", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,0.05), xlim=range(gen_sub[,8]))
	for (i in 1:20)
		{
		rep_sub <- subset(gen_sub, gen_sub[,1]==i)
		junc_dens <- density(rep_sub[,8])
		points(junc_dens$x, junc_dens$y, type="l", lwd=1.5, col="dodgerblue")
		polygon(junc_dens, col=adjustcolor("dodgerblue", alpha.f=0.3), border=FALSE)
		}

	## q densities
	plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,50), xlim=c(0,1))
	for (i in 1:20)
		{
		rep_sub <- subset(gen_sub, gen_sub[,1]==i)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=1.5, col="dodgerblue")
		polygon(q_dens, col=adjustcolor("dodgerblue", alpha.f=0.3), border=FALSE)
		}

	## Q densities
	plot(0, type="n", xlab="Q", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,30), xlim=c(0,1))
	mtext(paste0("trial: ", trial_in, "; m: ", m_in), line=-1.9, adj=0.98)
	mtext(paste0("c: ", c_in, "; gen: ", gen_in), line=-3.15, adj=0.98)
	for (i in 1:20)
		{
		rep_sub <- subset(gen_sub, gen_sub[,1]==i)
		Q_dens <- density(rep_sub[,6])
		points(Q_dens$x, Q_dens$y, type="l", lwd=1.5, col="dodgerblue")
		polygon(Q_dens, col=adjustcolor("dodgerblue", alpha.f=0.3), border=FALSE)
		}

	dev.off()
	}





########################################
## 10) pca
########################################


library("vegan")
uniq_runs <- c("dmi", "dmi_e", "path", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)

gen <- 100
quartz(height=12, width=18)
par(mar=c(5,5,0,0), mfrow=c(4,6), oma=c(0,0,4,4))
ctr <- 0
for (i in 1:length(uniq_runs))
	{
	for (j in 1:length(uniq_m))
		{
		for (k in 1:length(uniq_c))
			{
			ctr <- ctr + 1
			if (uniq_c[k]==0)	{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen", gen, ".csv") }
			else					{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen", gen, ".csv") }
			dat_in <- read.csv(handle, header=FALSE)
			pca_out <- prcomp(dat_in[,9:518], center=TRUE, scale=FALSE)
			pc1_exp <- round(summary(pca_out)$importance[2,1]*100, 2)
			pc2_exp <- round(summary(pca_out)$importance[2,2]*100, 2)
			plot(pca_out$x[,1], pca_out$x[,2], type="n", xlab=paste0("PC1 (", pc1_exp, "%)"), ylab=paste0("PC2 (", pc2_exp, "%)"), cex.lab=1.75, cex.axis=1.5, las=1, xlim=c(min(pca_out$x[,1]*2), max(pca_out$x[,1]*2)), ylim=c(min(pca_out$x[,2]*1.1), max(pca_out$x[,2]*1.1))); box(lwd=2)
			pca_plot_data <- cbind(dat_in[,1:4], pca_out$x[,1:2])
			points(pca_out$x[,1], pca_out$x[,2], col="gray", cex=0.5)
			ordiellipse(pca_out, pca_plot_data[,1], conf=0.9, col="red")
			if (ctr < 7)
				{
				mtext(paste0("m = ", uniq_m[j]), line=2, cex=1.5)
				mtext(paste0("c = ", uniq_c[k]), line=0.25, cex=1.5)
				}
			if (ctr %% 6 == 0)
				{
				if 		(uniq_runs[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
				}
			}
		}
	}



########################################
## 11) anova plots
########################################

hexagon <- function(x_center=NA, y_center=NA, scale=NA, plot_xmin=NA, plot_xmax=NA, plot_ymin=NA, plot_ymax=NA, color=NA)
	## function allows you to plot fillable hexagon points in plots
	## x_center and y_center mark the center of the point
	## scale controls the size of the point relative to the size of the plotting area (e.g., 0.1)
	## plot_xmin, plot_xmax, plot_ymin, and plot_ymax need to be the xlims and ylims of your plot
	## color is any color you'd like to use
	{
	x_scale <- scale * (plot_xmax - plot_xmin)
	y_scale <- scale * (plot_ymax - plot_ymin)
	ax <- x_center + x_scale
	ay <- y_center
	bx <- x_center + x_scale / 2
	by <- y_center + sqrt(3) * y_scale / 2
	cx <- x_center - x_scale / 2
	cy <- y_center + sqrt(3) * y_scale / 2
	dx <- x_center - x_scale
	dy <- y_center
	ex <- x_center - x_scale / 2
	ey <- y_center - sqrt(3) * y_scale / 2
	fx <- x_center + x_scale / 2
	fy <- y_center - sqrt(3) * y_scale / 2
	x_points <- c(ax, bx, cx, dx, ex, fx)
	y_points <- c(ay, by, cy, dy, ey, fy)
	polygon(x_points, y_points, col=color, border="black")
	}



uniq_runs <- c("dmi", "dmi_e", "path", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)
uniq_gen <- c(10, 100)

anova_out <- matrix(0, 48, 13)

ctr <- 1
for (i in 1:length(uniq_runs))
	{
	for (j in 1:length(uniq_m))
		{
		for (k in 1:length(uniq_c))
			{
			for (l in 1:length(uniq_gen))
				{
				print(ctr)
				if (uniq_c[k]==0)	{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen", uniq_gen[l], ".csv") }
				else							{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen", uniq_gen[l], ".csv") }
				dat_in <- read.csv(handle, header=FALSE)
				q_anova <- aov(dat_in[,5] ~ as.factor(dat_in[,1]))
				het_anova <- aov(dat_in[,6] ~ as.factor(dat_in[,1]))
				junc_anova <- aov(dat_in[,8] ~ as.factor(dat_in[,1]))
				
				anova_out[ctr,1] <- uniq_runs[i]
				anova_out[ctr,2] <- uniq_m[j]
				anova_out[ctr,3] <- uniq_c[k]
				anova_out[ctr,4] <- uniq_gen[l]
				anova_out[ctr,5] <- log10(summary(q_anova)[[1]][1,4])
				anova_out[ctr,6] <- log10(summary(het_anova)[[1]][1,4])
				anova_out[ctr,7] <- log10(summary(junc_anova)[[1]][1,4])
				anova_out[ctr,8] <- summary(q_anova)[[1]][1,4]
				anova_out[ctr,9] <- summary(q_anova)[[1]][1,5]
				anova_out[ctr,10] <- summary(het_anova)[[1]][1,4]
				anova_out[ctr,11] <- summary(het_anova)[[1]][1,5]
				anova_out[ctr,12] <- summary(junc_anova)[[1]][1,4]
				anova_out[ctr,13] <- summary(junc_anova)[[1]][1,5]
				ctr <- ctr+1
				}
			}
		}
	}
	
## export for table
table_anova_out <- anova_out[,c(1:4,8:13)]
colnames(table_anova_out) <- c("run", "m", "c", "gen", "q_f", "q_p", "het_f", "het_p", "junc_f", "junc_p")
#write.table(table_anova_out, file="anova_out_table.txt", row.names=F, quote=F, sep="\t")

anova_out_gen10 <- subset(anova_out, anova_out[,4]=="10")
anova_out_gen100 <- subset(anova_out, anova_out[,4]=="100")


anova_cols <- met.brewer("Hokusai1", 24)


## correlations

## gen10
cor(as.numeric(anova_out_gen10[,5]), as.numeric(anova_out_gen10[,6]))
	## [1] 0.8739796

## gen100
cor(as.numeric(anova_out_gen100[,5]), as.numeric(anova_out_gen100[,6]))
	## [1] 0.9332555




## gen 10 side by side

quartz(height=8, width=20)
layout(matrix(c(25,25,25,25,1,2,3,4,5,6,
				25,25,25,25,7,8,9,10,11,12,
				25,25,25,25,13,14,15,16,17,18,
				25,25,25,25,19,20,21,22,23,24), 4, 10, byrow=TRUE))
par(mar=c(5,5,0,0), oma=c(0,0,4,4))

for (i in 1:24)
	{
	if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	else										{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	dat_in <- read.csv(handle, header=FALSE)
	plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1)
	if		(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.2)		{ points(0.9, 10.8, pch=21, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.9)		{ points(0.9, 10.8, pch=22, bg=anova_cols[i], cex=3) }
	#else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]=="neutral")	{ points(0.9, 10.8, pch=13, col=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0)	{ hexagon(x_center=0.9, y_center=10.8, scale=0.075, plot_xmin=0, plot_xmax=1, plot_ymin=0, plot_ymax=12, color=anova_cols[i]) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.2)			{ points(0.9, 10.8, pch=24, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.9)			{ points(0.9, 10.8, pch=25, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0)	{ points(0.9, 10.8, pch=23, bg=anova_cols[i], cex=3) }
	for (j in 1:20)
		{
		rep_sub <- subset(dat_in, dat_in[,1]==j)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=0.5, col=anova_cols[i])
		polygon(q_dens, col=adjustcolor(anova_cols[i], alpha.f=0.2), border=FALSE)
		}
	if (i < 7)
		{
		mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
		mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
		}
	if (i %% 6 == 0)
		{
		if (anova_out_gen10[i,1]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
		}
	}
plot(anova_out_gen10[,5], anova_out_gen10[,6], type="n", xlab="q anova log10 F", ylab="Q anova log10 F", cex.lab=2.5, cex.axis=1.5, las=1, xlim=range(as.numeric(anova_out_gen10[,5:6])), ylim=range(as.numeric(anova_out_gen10[,5:6])))
abline(0, 1, lty=2, lwd=4)
for (i in 1:24)
	{
	if		(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.2)		{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=21, bg=anova_cols[i], cex=5) }
	else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.9)		{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=22, bg=anova_cols[i], cex=4) }
	else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0)	{ hexagon(as.numeric(anova_out_gen10[i,5]), as.numeric(anova_out_gen10[i,6]), scale=0.02, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color=anova_cols[i]) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.2)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=24, bg=anova_cols[i], cex=4) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.9)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=25, bg=anova_cols[i], cex=4) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0)	{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=23, bg=anova_cols[i], cex=4) }
	}
legend("topleft", legend=c("m = 0.01; c = 0", "m = 0.01; c = 0.2", "m = 0.01; c = 0.9", "m = 0.2; c = 0", "m = 0.2; c = 0.2", "m = 0.2; c = 0.9"), pch=c(20,21,22,23,24,25), pt.cex=3, pt.bg="black", cex=2)
hexagon(-0.178, 3.06, scale=0.015, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color="black")






## gen 100 side by side

quartz(height=8, width=20)
layout(matrix(c(25,25,25,25,1,2,3,4,5,6,
				25,25,25,25,7,8,9,10,11,12,
				25,25,25,25,13,14,15,16,17,18,
				25,25,25,25,19,20,21,22,23,24), 4, 10, byrow=TRUE))
par(mar=c(5,5,0,0), oma=c(0,0,4,4))

for (i in 1:24)
	{
	if (anova_out_gen100[i,3]==0)	{ handle <- paste0(anova_out_gen100[i,1], "_m", anova_out_gen100[i,2], "_neutral", "_deme6_gen", anova_out_gen100[i,4], ".csv") }
	else										{ handle <- paste0(anova_out_gen100[i,1], "_m", anova_out_gen100[i,2], "_c", anova_out_gen100[i,3], "_deme6_gen", anova_out_gen100[i,4], ".csv") }
	dat_in <- read.csv(handle, header=FALSE)
	
	plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,50), xlim=c(0,1), las=1)
	if		(anova_out_gen100[i,2]==0.01 & anova_out_gen100[i,3]==0.2)			{ points(0.9, 45, pch=21, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen100[i,2]==0.01 & anova_out_gen100[i,3]==0.9)			{ points(0.9, 45, pch=22, bg=anova_cols[i], cex=3) }
	#else if	(anova_out_gen100[i,2]==0.01 & anova_out_gen100[i,3]=="neutral")	{ points(0.9, 45, pch=13, col=anova_cols[i], cex=3) }
	else if	(anova_out_gen100[i,2]==0.01 & anova_out_gen10[i,3]==0)		{ hexagon(x_center=0.9, y_center=45, scale=0.075, plot_xmin=0, plot_xmax=1, plot_ymin=0, plot_ymax=50, color=anova_cols[i]) }
	else if	(anova_out_gen100[i,2]==0.2 & anova_out_gen100[i,3]==0.2)			{ points(0.9, 45, pch=24, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen100[i,2]==0.2 & anova_out_gen100[i,3]==0.9)			{ points(0.9, 45, pch=25, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen100[i,2]==0.2 & anova_out_gen100[i,3]==0)		{ points(0.9, 45, pch=23, bg=anova_cols[i], cex=3) }
	for (j in 1:20)
		{
		rep_sub <- subset(dat_in, dat_in[,1]==j)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=0.5, col=anova_cols[i])
		polygon(q_dens, col=adjustcolor(anova_cols[i], alpha.f=0.2), border=FALSE)
		}
	if (i < 7)
		{
		mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
		mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
		}
	if (i %% 6 == 0)
		{
		if (anova_out_gen10[i,1]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
		}
	}
	
plot(anova_out_gen100[,5], anova_out_gen100[,6], type="n", xlab="q anova log10 F", ylab="Q anova log10 F", cex.lab=2.5, cex.axis=1.5, las=1, xlim=range(as.numeric(anova_out_gen100[,5:6])), ylim=range(as.numeric(anova_out_gen100[,5:6])))
abline(0, 1, lty=2, lwd=4)
for (i in 1:24)
	{
	if		(anova_out_gen100[i,2]==0.01 & anova_out_gen100[i,3]==0.2)			{ points(anova_out_gen100[i,5], anova_out_gen100[i,6], pch=21, bg=anova_cols[i], cex=5) }
	else if	(anova_out_gen100[i,2]==0.01 & anova_out_gen100[i,3]==0.9)			{ points(anova_out_gen100[i,5], anova_out_gen100[i,6], pch=22, bg=anova_cols[i], cex=4) }
	else if	(anova_out_gen100[i,2]==0.01 & anova_out_gen10[i,3]==0)	{ hexagon(as.numeric(anova_out_gen100[i,5]), as.numeric(anova_out_gen100[i,6]), scale=0.02, plot_xmin=min(as.numeric(anova_out_gen100[,5:6])), plot_xmax=max(as.numeric(anova_out_gen100[,5:6])), plot_ymin=min(as.numeric(anova_out_gen100[,5:6])), plot_ymax=max(as.numeric(anova_out_gen100[,5:6])), color=anova_cols[i]) }
	else if	(anova_out_gen100[i,2]==0.2 & anova_out_gen100[i,3]==0.2)			{ points(anova_out_gen100[i,5], anova_out_gen100[i,6], pch=24, bg=anova_cols[i], cex=4) }
	else if	(anova_out_gen100[i,2]==0.2 & anova_out_gen100[i,3]==0.9)			{ points(anova_out_gen100[i,5], anova_out_gen100[i,6], pch=25, bg=anova_cols[i], cex=4) }
	else if	(anova_out_gen100[i,2]==0.2 & anova_out_gen100[i,3]==0)		{ points(anova_out_gen100[i,5], anova_out_gen100[i,6], pch=23, bg=anova_cols[i], cex=4) }
	}
legend("topleft", legend=c("m = 0.01; c = 0", "m = 0.01; c = 0.2", "m = 0.01; c = 0.9", "m = 0.2; c = 0", "m = 0.2; c = 0.2", "m = 0.2; c = 0.9"), pch=c(20,21,22,23,24,25), pt.cex=3, pt.bg="black", cex=2)
hexagon(-0.06, 3.49, scale=0.015, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color="black")
box(lwd=2)



## gen10 right panels for paper (different colors)

right_panel_cols <- met.brewer("Kandinsky", 3)
right_panel_cols_list <- rep(rev(right_panel_cols), 8)

quartz(height=8, width=12)
par(mar=c(5,5,0,0), oma=c(0,0,4,4), mfrow=c(4,6))

for (i in 1:24)
	{
	if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	else							{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	dat_in <- read.csv(handle, header=FALSE)
	plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.75, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1, xaxt="n")
	axis(1, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=1.5)
	for (j in 1:20)
		{
		rep_sub <- subset(dat_in, dat_in[,1]==j)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=0.5, col=right_panel_cols_list[i])
		polygon(q_dens, col=adjustcolor(right_panel_cols_list[i], alpha.f=0.2), border=FALSE)
		}
	if (i < 7)
		{
		mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
		mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
		}
	if (i %% 6 == 0)
		{
		if (anova_out_gen10[i,1]=="dmi") { mtext("DMI", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="dmi_e") { mtext("DMI + env", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path") { mtext("Path", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path_e") { mtext("Path + env", side=4, line=1.25, cex=1.5) }
		}
	}



## test different colors for different replicates

triangle_cols <- met.brewer("OKeeffe1", 20)

quartz(height=8, width=12)
par(mar=c(5,5,0,0), oma=c(0,0,4,4), mfrow=c(4,6))

for (i in 1:24)
	{
	if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	else							{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	dat_in <- read.csv(handle, header=FALSE)
	plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.75, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1, xaxt="n")
	rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col="light gray")
	axis(1, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=1.5)
	for (j in 1:20)
		{
		rep_sub <- subset(dat_in, dat_in[,1]==j)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=1.5, col= triangle_cols[j])
		}
	if (i < 7)
		{
		mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
		mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
		}
	if (i %% 6 == 0)
		{
		if (anova_out_gen10[i,1]=="dmi") { mtext("DMI", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="dmi_e") { mtext("DMI + env", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path") { mtext("Path", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path_e") { mtext("Path + env", side=4, line=1.25, cex=1.5) }
		}
	}






## gen10 talk right panels (full)

quartz(height=8, width=12)
par(mar=c(5,5,0,0), oma=c(0,0,4,4), mfrow=c(4,6))

for (i in 1:24)
	{
	if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	else										{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme6_gen", anova_out_gen10[i,4], ".csv") }
	dat_in <- read.csv(handle, header=FALSE)
	plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1)
	for (j in 1:20)
		{
		rep_sub <- subset(dat_in, dat_in[,1]==j)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=0.5, col=anova_cols[i])
		polygon(q_dens, col=adjustcolor(anova_cols[i], alpha.f=0.2), border=FALSE)
		}
	if (i < 7)
		{
		mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
		mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
		}
	if (i %% 6 == 0)
		{
		if (anova_out_gen10[i,1]=="dmi") { mtext("DMI", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="dmi_e") { mtext("DMI + env", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path") { mtext("Path", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path_e") { mtext("Path + env", side=4, line=1.25, cex=1.5) }
		}
	}


## gen10 talk left panel only

quartz(height=8, width=8)
par(mar=c(5,5,1,1))
plot(anova_out_gen10[,5], anova_out_gen10[,6], type="n", xlab="q anova log10 F", ylab="Q anova log10 F", cex.lab=2.5, cex.axis=1.5, las=1, xlim=range(as.numeric(anova_out_gen10[,5:6])), ylim=range(as.numeric(anova_out_gen10[,5:6])))
abline(0, 1, lty=2, lwd=4)
for (i in 1:24)
	{
	if		(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.2)		{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=21, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.9)		{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=22, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0)	{ hexagon(as.numeric(anova_out_gen10[i,5]), as.numeric(anova_out_gen10[i,6]), scale=0.02, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color=anova_cols[i]) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.2)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=24, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.9)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=25, bg=anova_cols[i], cex=3) }
	else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=23, bg=anova_cols[i], cex=3) }
	}
legend("topleft", legend=c("m = 0.01; c = 0", "m = 0.01; c = 0.2", "m = 0.01; c = 0.9", "m = 0.2; c = 0", "m = 0.2; c = 0.2", "m = 0.2; c = 0.9"), pch=c(20,21,22,23,24,25), pt.cex=2.5, pt.bg="black", cex=1.5)
hexagon(-0.16, 3.025, scale=0.02, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color="black")
box(lwd=2)






## example densities for talk

quartz(height=6, width=9)
par(mar=c(1,1,1,1), mfrow=c(2,3))
plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), type="n", xaxt="n", yaxt="n"); box(lwd=3)
	polygon(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), col="gray", border=NA)
	lines(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), type="l", lwd=7)
plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), type="n", xaxt="n", yaxt="n"); box(lwd=3)
	polygon(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), col="gray", border=NA)
	lines(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), type="l", lwd=7)
plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), type="n", xaxt="n", yaxt="n"); box(lwd=3)
	polygon(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), col="gray", border=NA)
	lines(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 20, 4)+dbeta(seq(0, 1, 0.01), 3, 10), type="l", lwd=7)

plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 5, 5)+dbeta(seq(0, 1, 0.01), 9, 30), type="n", xaxt="n", yaxt="n"); box(lwd=3)
	polygon(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 5, 5)+dbeta(seq(0, 1, 0.01), 9, 30), col="gray", border=NA)
	lines(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 5, 5)+dbeta(seq(0, 1, 0.01), 9, 30), type="l", lwd=7)
plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 40, 40), type="n", xaxt="n", yaxt="n"); box(lwd=3)
	polygon(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 40, 40), col="gray", border=NA)
	lines(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 40, 40), type="l", lwd=7)
plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 15, 3), type="n", xaxt="n", yaxt="n"); box(lwd=3)
	polygon(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 15, 3), col="gray", border=NA)
	lines(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), 15, 3), type="l", lwd=7)









########################################
## 12) contingency
########################################

## NOTE: dmis found at locations 1.4 and 1.6 (columns 12 & 14)


library(MetBrewer)
contin_cols <- met.brewer("Derain", 151)

contin <- read.csv("dmi_m0.2_c0.9_deme6_gen10.csv", header=F)
	dim(contin)
	contin[1:10,1:10]


## for 1.4 vs 1.6, columns should be set to 12 and 14
## for 2.4 vs 2.6, columns should be set to 63 and 65
## for 3.4 vs 3.6, columns should be set to 114 and 116

contin_mat <- matrix(0, 20, 9)
for (i in 1:20)
	{
	rep_sub <- subset(contin, contin[,1]==i)
	for (j in 1:dim(rep_sub)[1])
		{
		if (rep_sub[j,12]==0)
			{
			if		(rep_sub[j,14]==0) { contin_mat[i,1] <- contin_mat[i,1] + 1 }
			else if	(rep_sub[j,14]==1) { contin_mat[i,2] <- contin_mat[i,2] + 1 }
			else if	(rep_sub[j,14]==2) { contin_mat[i,3] <- contin_mat[i,3] + 1 }
			}
		else if (rep_sub[j,12]==1)
			{
			if		(rep_sub[j,14]==0) { contin_mat[i,4] <- contin_mat[i,4] + 1 }
			else if	(rep_sub[j,14]==1) { contin_mat[i,5] <- contin_mat[i,5] + 1 }
			else if	(rep_sub[j,14]==2) { contin_mat[i,6] <- contin_mat[i,6] + 1 }
			}
		if (rep_sub[j,12]==2)
			{
			if		(rep_sub[j,14]==0) { contin_mat[i,7] <- contin_mat[i,7] + 1 }
			else if	(rep_sub[j,14]==1) { contin_mat[i,8] <- contin_mat[i,8] + 1 }
			else if	(rep_sub[j,14]==2) { contin_mat[i,9] <- contin_mat[i,9] + 1 }
			}
		}
	}


quartz(height=6, width=14)
par(mar=c(1,3,3,1), mfrow=c(3,7))
for (i in 1:20)
	{
	plot(0, type="n", xlim=c(0,3), ylim=c(0,3), xaxt="n", yaxt="n", xlab="", ylab="")
	axis(2, at=c(0.5, 1.5, 2.5), labels=c(2, 1, 0), las=1, cex.axis=2)
	axis(3, at=c(0.5, 1.5, 2.5), labels=c(0, 1, 2), cex.axis=2)
	rect(0, 2, 1, 3, col=contin_cols[(contin_mat[i,1] + 1)], lwd=2)
	rect(1, 2, 2, 3, col=contin_cols[(contin_mat[i,2] + 1)], lwd=2)
	rect(2, 2, 3, 3, col=contin_cols[(contin_mat[i,3] + 1)], lwd=2)
	rect(0, 1, 1, 2, col=contin_cols[(contin_mat[i,4] + 1)], lwd=2)
	rect(1, 1, 2, 2, col=contin_cols[(contin_mat[i,5] + 1)], lwd=2)
	rect(2, 1, 3, 2, col=contin_cols[(contin_mat[i,6] + 1)], lwd=2)
	rect(0, 0, 1, 1, col=contin_cols[(contin_mat[i,7] + 1)], lwd=2)
	rect(1, 0, 2, 1, col=contin_cols[(contin_mat[i,8] + 1)], lwd=2)
	rect(2, 0, 3, 1, col=contin_cols[(contin_mat[i,9] + 1)], lwd=2)
	}

legend_image <- as.raster(matrix(rev(contin_cols), ncol=1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
rasterImage(legend_image, 1/3, 0.05, 2/3, 0.95)
rect(1/3, 0.05, 2/3, 0.95, lwd=2)
text(0.85, 0.05, labels="0", cex=2)
text(0.85, 0.5, labels="75", cex=2)
text(0.85, 0.95, labels="150", cex=2)



########################################
## 13) qQ across replicates
########################################

library(MetBrewer)
triangle_cols <- met.brewer("OKeeffe1", 20)

triangle <- read.csv("dmi_m0.01_c0.9.main_first8cols.txt", header=TRUE)
	dim(triangle)
	head(triangle)

gen10 <- subset(triangle, triangle[,2]==10)
	dim(gen10)
	head(gen10)

gen10_deme6 <- subset(gen10, gen10[,3]==6)
	dim(gen10_deme6)
	head(gen10_deme6)

quartz(height=12, width=15)
par(mar=c(5,5,1,1), mfrow=c(4,5))
for (i in 1:20)
	{
	rep_sub <- subset(gen10_deme6, gen10_deme6[,1]==i)
	plot(0, type="n", xlab="q", ylab="Q", xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
	segments(0, 0, 0.5, 1, lwd=3)
	segments(1, 0, 0.5, 1, lwd=3)
	for (j in 1:150)
		{
		points(rep_sub[j,5], rep_sub[j,6], pch=21, bg=adjustcolor(triangle_cols[i], alpha.f=0.75), cex=2)
		}
	box(lwd=2)
	}





########################################
## 14) gen10 vs gen100
########################################

library(MetBrewer)
rep_cols <- met.brewer("OKeeffe1", 20)

## q

uniq_runs <- c("dmi", "dmi_e", "path", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)

quartz(height=12, width=18)
par(mar=c(5,5,1,1), mfrow=c(4,6), oma=c(0,0,4,4))
ctr <- 0
for (i in 1:length(uniq_runs))
	{
	for (j in 1:length(uniq_m))
		{
		for (k in 1:length(uniq_c))
			{
			ctr <- ctr + 1
			if (uniq_c[k]==0)
				{
				handle10 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen10.csv")
				handle100 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen100.csv")
				}
			else	
				{
				handle10 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen10.csv")
				handle100 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen100.csv")
				}
				
			dat10 <- read.csv(handle10, header=FALSE)
			dat100 <- read.csv(handle100, header=FALSE)
			plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="Median gen 10 q", ylab="Median gen 100 q", cex.lab=1.75, cex.axis=1.5, las=1)
			box(lwd=2)
			if (ctr < 7)
				{
				mtext(paste0("m = ", uniq_m[j]), line=2, cex=1.5)
				mtext(paste0("c = ", uniq_c[k]), line=0.25, cex=1.5)
				}
			if (ctr %% 6 == 0)
				{
				if 		(uniq_runs[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
				}
			cor_mat <- matrix(0, 20, 2)
			for (l in 1:20)
				{
				rep_sub10 <- subset(dat10, dat10[,1]==l)
				rep_sub100 <- subset(dat100, dat100[,1]==l)
				cor_mat[l,1] <- median(rep_sub10[,5])
				cor_mat[l,2] <- median(rep_sub100[,5])
				points(median(rep_sub10[,5]), median(rep_sub100[,5]), pch=21, bg=rep_cols[l], cex=2)
				}
			mtext(paste0("r = ", round(cor(cor_mat[,1], cor_mat[,2], method="pearson"), 3)), line=-2, adj=0.05, cex=1.25)
			}
		}
	}




## Q

uniq_runs <- c("dmi", "dmi_e", "path", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)

quartz(height=12, width=18)
par(mar=c(5,5,1,1), mfrow=c(4,6), oma=c(0,0,4,4))
ctr <- 0
for (i in 1:length(uniq_runs))
	{
	for (j in 1:length(uniq_m))
		{
		for (k in 1:length(uniq_c))
			{
			ctr <- ctr + 1
			if (uniq_c[k]==0)
				{
				handle10 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen10.csv")
				handle100 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen100.csv")
				}
			else	
				{
				handle10 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen10.csv")
				handle100 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen100.csv")
				}
				
			dat10 <- read.csv(handle10, header=FALSE)
			dat100 <- read.csv(handle100, header=FALSE)
			plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="Median gen 10 Q", ylab="Median gen 100 Q", cex.lab=1.75, cex.axis=1.5, las=1)
			box(lwd=2)
			if (ctr < 7)
				{
				mtext(paste0("m = ", uniq_m[j]), line=2, cex=1.5)
				mtext(paste0("c = ", uniq_c[k]), line=0.25, cex=1.5)
				}
			if (ctr %% 6 == 0)
				{
				if 		(uniq_runs[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
				else if (uniq_runs[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
				}
			cor_mat <- matrix(0, 20, 2)
			for (l in 1:20)
				{
				rep_sub10 <- subset(dat10, dat10[,1]==l)
				rep_sub100 <- subset(dat100, dat100[,1]==l)
				cor_mat[l,1] <- median(rep_sub10[,6])
				cor_mat[l,2] <- median(rep_sub100[,6])
				points(median(rep_sub10[,6]), median(rep_sub100[,6]), pch=21, bg=rep_cols[l], cex=2)
				}
			mtext(paste0("r = ", round(cor(cor_mat[,1], cor_mat[,2], method="pearson"), 3)), line=-2, adj=0.05, cex=1.25)
			}
		}
	}






















