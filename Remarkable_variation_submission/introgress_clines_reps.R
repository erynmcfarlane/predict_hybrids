

##### genomic clines across replicates
source("SNP_inputs.R")

source("genomic.cline.plot.reps.R")

###I'm not sure it makes sense using all 510 snps for this. let's do just the 3?
data_deme_6_gen_10[,c(1:8, 519:521, 12, 18, 114)]->alldata_df_6_10
alldata_df_6_10$mech<-as.factor(ifelse(alldata_df_6_10$mech=="dmi_m", "dmi", ifelse(alldata_df_6_10$mech=="path_m", "path",ifelse(alldata_df_6_10$mech=="path_e", 'path_e', "dmi_e"))))

alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "path_e")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "path")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "dmi_e")
alldata_df_6_10$mech<-relevel(alldata_df_6_10$mech, "dmi")

### I want to do this separately for the 24 categories we have! ###
alldata_df_6_10$index<-paste(alldata_df_6_10$m, alldata_df_6_10$c, alldata_df_6_10$mech)
alldata_df_6_10$index_reps<-paste(alldata_df_6_10$m, alldata_df_6_10$c, alldata_df_6_10$mech, alldata_df_6_10$rep)
alldata_df_6_10$index<-as.factor(alldata_df_6_10$index)

alldata_df_6_10_noE<-alldata_df_6_10[which(alldata_df_6_10$mech %in% c("dmi", "path")),]
alldata_df_6_10_noE$rep<-as.factor(alldata_df_6_10_noE$rep)
colours<-met.brewer(name='OKeeffe1', n=20, type='continuous') 

### Figure S8###
pdf(file="genomic_cline_plots1.4.pdf", width=25, height=10)
par(mfrow=c(2,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(9,7,8,12,10,11,
                3,1,2,6,4,5), 2, 6, byrow=TRUE))
for(i in 1:length(unique(alldata_df_6_10_noE$index))){
  genomic.clines.reps<-list()
  Fitted.AA<-list()
  Fitted.Aa<-list()
    for(j in 1:20){
  introgress.data<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),c(12:14)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),]$q
  loci.data<-matrix(nrow=length(introgress.data), ncol=3)
  dim(loci.data)<-c(3, 3)
  loci.data[,1]<-colnames(introgress.data)
  loci.data[,2]<-'C'
  loci.data[,3]<-chromosome
  genomic.clines.reps[[j]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=hi.index, loci.data=loci.data)
  Fitted.AA[[j]]<-t(genomic.clines.reps[[j]]$Fitted.AA)
  Fitted.Aa[[j]]<-t(genomic.clines.reps[[j]]$Fitted.Aa)
    }
  Fitted.AA.df<-do.call(rbind.data.frame, Fitted.AA)
  Fitted.AA.df$rep<-as.factor(rep(seq(1,20), each=150))
  Fitted.Aa.df<-do.call(rbind.data.frame, Fitted.Aa)
  Fitted.Aa.df$rep<-as.factor(rep(seq(1,20), each=150))
  

  
  plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pty='s', xaxt="n", yaxt="n")
  axis(1, at=c(0,1), cex.axis=2.5)
  axis(2, at=c(0,1), cex.axis=2.5, las=1)
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
  genomic.cline.plot(genomic.clines.reps, 1)
  if (i %in% c(9,7,8,12,10,11))
    {
    mtext(paste0("m = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$m[i]), line=2, cex=1.5)
    mtext(paste0("c = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$c[i]), line=0.25, cex=1.5)
    }
  
  if (i %in% c(11,5))
  {
    if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
    }
}

mtext('Admixture Proportion', side = 1, outer = TRUE, line = 2, cex=2.5)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2.5)
dev.off()

####Figure S9 ####

pdf(file="genomic_cline_plots1.10.pdf", width=25, height=10)
par(mfrow=c(2,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(9,7,8,12,10,11,
                3,1,2,6,4,5), 2, 6, byrow=TRUE))
for(i in 1:length(unique(alldata_df_6_10_noE$index))){
  genomic.clines.reps<-list()
  Fitted.AA<-list()
  Fitted.Aa<-list()
  for(j in 1:20){
    introgress.data<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),c(12:14)]
    chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
    hi.index<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),]$q
    loci.data<-matrix(nrow=length(introgress.data), ncol=3)
    dim(loci.data)<-c(3, 3)
    loci.data[,1]<-colnames(introgress.data)
    loci.data[,2]<-'C'
    loci.data[,3]<-chromosome
    genomic.clines.reps[[j]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=hi.index, loci.data=loci.data)
    Fitted.AA[[j]]<-t(genomic.clines.reps[[j]]$Fitted.AA)
    Fitted.Aa[[j]]<-t(genomic.clines.reps[[j]]$Fitted.Aa)
  }
  Fitted.AA.df<-do.call(rbind.data.frame, Fitted.AA)
  Fitted.AA.df$rep<-as.factor(rep(seq(1,20), each=150))
  Fitted.Aa.df<-do.call(rbind.data.frame, Fitted.Aa)
  Fitted.Aa.df$rep<-as.factor(rep(seq(1,20), each=150))
  
  
  
  plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pty='s', xaxt="n", yaxt="n")
  axis(1, at=c(0,1), cex.axis=2.5)
  axis(2, at=c(0,1), cex.axis=2.5, las=1)
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
  genomic.cline.plot(genomic.clines.reps, 2) ###this is where I tell it which SNP of the three I want
  if (i %in% c(9,7,8,12,10,11))
  {
    mtext(paste0("m = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$m[i]), line=2, cex=1.5)
    mtext(paste0("c = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$c[i]), line=0.25, cex=1.5)
  }
  
  if (i %in% c(11,5))
  {
    if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
  }
}

mtext('Admixture Proportion', side = 1, outer = TRUE, line = 2, cex=2.5)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2.5)
dev.off()

###Figure S10

pdf(file="genomic_cline_plots3.4.pdf", width=25, height=10)
par(mfrow=c(2,6), mar=c(5,5,0,0), oma=c(5,5,4,4))
layout(matrix(c(9,7,8,12,10,11,
                3,1,2,6,4,5), 2, 6, byrow=TRUE))
for(i in 1:length(unique(alldata_df_6_10_noE$index))){
  genomic.clines.reps<-list()
  Fitted.AA<-list()
  Fitted.Aa<-list()
  for(j in 1:20){
    introgress.data<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),c(12:14)]
    chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
    hi.index<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i] & alldata_df_6_10_noE$rep==j),]$q
    loci.data<-matrix(nrow=length(introgress.data), ncol=3)
    dim(loci.data)<-c(3, 3)
    loci.data[,1]<-colnames(introgress.data)
    loci.data[,2]<-'C'
    loci.data[,3]<-chromosome
    genomic.clines.reps[[j]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=hi.index, loci.data=loci.data)
    Fitted.AA[[j]]<-t(genomic.clines.reps[[j]]$Fitted.AA)
    Fitted.Aa[[j]]<-t(genomic.clines.reps[[j]]$Fitted.Aa)
  }
  Fitted.AA.df<-do.call(rbind.data.frame, Fitted.AA)
  Fitted.AA.df$rep<-as.factor(rep(seq(1,20), each=150))
  Fitted.Aa.df<-do.call(rbind.data.frame, Fitted.Aa)
  Fitted.Aa.df$rep<-as.factor(rep(seq(1,20), each=150))
  
  
  
  plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pty='s', xaxt="n", yaxt="n")
  axis(1, at=c(0,1), cex.axis=2.5)
  axis(2, at=c(0,1), cex.axis=2.5, las=1)
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
  genomic.cline.plot(genomic.clines.reps, 3) ###this is where I tell it which SNP of the three I want
  if (i %in% c(9,7,8,12,10,11))
  {
    mtext(paste0("m = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$m[i]), line=2, cex=1.5)
    mtext(paste0("c = ", alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$c[i]), line=0.25, cex=1.5)
  }
  
  if (i %in% c(11,5))
  {
    if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
    else if (alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),]$mech[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
  }
}

mtext('Admixture Proportion', side = 1, outer = TRUE, line = 2, cex=2.5)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 2, cex=2.5)
dev.off()



####Figure 4 ####
pdf(file="genomic_cline_plots_0.01_0.09_DMI.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4), mpg=c(3,3,0))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot - ERYN, rewrite the function so you can put loci number in it at some point.
for(j in 1:20){
  introgress.data<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index=="0.01 0.9 dmi" & alldata_df_6_10_noE$rep==j),c(12:14)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_6_10_noE[which(alldata_df_6_10_noE$index=="0.01 0.9 dmi" & alldata_df_6_10_noE$rep==j),]$q
  loci.data<-matrix(nrow=length(introgress.data), ncol=3)
  dim(loci.data)<-c(3, 3)
  loci.data[,1]<-colnames(introgress.data)
  loci.data[,2]<-'C'
  loci.data[,3]<-chromosome
  genomic.clines.reps[[j]]<-genomic.clines(introgress.data=t(introgress.data), hi.index=hi.index, loci.data=loci.data)
  Fitted.AA[[j]]<-t(genomic.clines.reps[[j]]$Fitted.AA)
  Fitted.Aa[[j]]<-t(genomic.clines.reps[[j]]$Fitted.Aa)
}
##SNP 1.4
plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pty='s', xaxt="n", yaxt="n", mar=c(1,1,0,0))
axis(1, at=c(0,1), cex.axis=3)
axis(2, at=c(0,1), cex.axis=3, las=1)
rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
genomic.cline.plot(genomic.clines.reps, 1)
##SNP 1.10
plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pty='s', xaxt="n", yaxt="n")
axis(1, at=c(0,1), cex.axis=3)
axis(2, at=c(0,1), cex.axis=3, las=1)
rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
genomic.cline.plot(genomic.clines.reps, 2)
##SNP 3.4
plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pty='s', xaxt="n", yaxt="n")
axis(1, at=c(0,1), cex.axis=3)
axis(2, at=c(0,1), cex.axis=3, las=1)
rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col='light gray')
genomic.cline.plot(genomic.clines.reps, 3)

mtext("Selected SNP", line=1, cex=4, at=-1.85)
mtext("LD with selected SNP", line=1, cex=4, at=-0.68)
mtext("Neutral SNP", line=1, cex=4, at=0.5)
mtext('Admixture Proportion', side = 1, outer = TRUE, line = 1.5, cex=4)
mtext('Probability of Genotype', side = 2, outer = TRUE, line = 1, cex=4)

dev.off()

### Table S4
l1.4_genotype_fstat<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l1.4_genotype_pvalue<-vector(length = length(unique(alldata_df_6_10_noE$index)))

l1.10_genotype_fstat<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l1.10_genotype_pvalue<-vector(length = length(unique(alldata_df_6_10_noE$index)))

l3.4_genotype_fstat<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l3.4_genotype_pvalue<-vector(length = length(unique(alldata_df_6_10_noE$index)))


for(i in 1:length(unique(alldata_df_6_10_noE$index))){
l1.4_genotype_fstat[[i]]<-unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[7]
l1.4_genotype_pvalue[[i]]<-unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[9]

l1.10_genotype_fstat[[i]]<-unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[7]
l1.10_genotype_pvalue[[i]]<-unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[9]

l3.4_genotype_fstat[[i]]<-unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[7]
l3.4_genotype_pvalue[[i]]<-unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[9]
} 

Genotype.Anova<-cbind(unique(as.character(alldata_df_6_10_noE$index)), l1.4_genotype_fstat, l1.4_genotype_pvalue, l1.10_genotype_fstat, l1.10_genotype_pvalue, l3.4_genotype_fstat, l3.4_genotype_pvalue)  
 write.table(Genotype.Anova, file="Genotype.Anova.repfactor.csv", col.names=T, row.names=F, quote=F, sep=',') 
  