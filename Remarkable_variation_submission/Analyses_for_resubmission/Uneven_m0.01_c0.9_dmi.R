#### between here and the next section is what's needed to run the asymetrical simulations
module load arcc/1.0
module load gcc/12.2.0
module load gsl/2.7.1
module load r

/gscratch/buerkle/data/incompatible/dfuse_src/dfuse


/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0.9 -o dmi_m0.01_0.9_uneven

#### once all of the simulations are done, move into R to do all this ####
library(tidyverse)
library(MASS)
library(data.table)
library(MetBrewer)
library(introgress)

datafiles<-list.files("/gscratch/emcfarl2/predicting_hybrids/",  pattern="*uneven.main", recursive=TRUE, include.dirs=TRUE)
basenames<-basename(datafiles)
m<-str_extract(basenames, "(\\d+\\.*\\d*)")
c<-str_match(basenames,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(basenames, "^([^_]+_){1}([^_])") 

alldata<-list()
for(i in 1:length(datafiles)){
  alldata[[i]]<-fread(datafiles[i], sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}


alldata_df<-do.call(rbind.data.frame, alldata)

rm(alldata)

alldata_df[which(alldata_df$deme==2), ]->data_deme_2
data_deme_2[which(data_deme_2$gen==10),]->data_deme_2_gen_10


genomic.cline.plot<-function(cline.data, snp.no){ ### I've now put this in to tell it which SNP i want it to plot so that I don't have to have it hard coded in.
  for(k in 1:length(cline.data)){
    n.loci <- dim(cline.data[[k]]$Loci.data)[1]
    hi <- cline.data[[k]]$hybrid.index
    col<-colours 
    
    for (l in snp.no:snp.no) {
      AA.line <- cline.data[[k]]$Fitted.AA[l, ]
      line.matrix <- rbind(hi, AA.line)[, order(hi)]
      lines(line.matrix[1, ], line.matrix[2, ], lty = 1,lwd=4, col=col[k])
      
      Aa.line <- cline.data[[k]]$Fitted.Aa[l, ]
      line.matrix <- rbind(hi, Aa.line)[, order(hi)]
      lines(line.matrix[1, ], line.matrix[2, ], lty = 2, lwd=4,col=col[k])
    }
  }
}

data_deme_2_gen_10[,c(1:8, 519:521, 12, 18, 114)]->alldata_df_2_10
alldata_df_2_10$mech<-as.factor(ifelse(alldata_df_2_10$mech=="dmi_m", "dmi", ifelse(alldata_df_2_10$mech=="path_m", "path",ifelse(alldata_df_2_10$mech=="path_e", 'path_e', "dmi_e"))))

alldata_df_2_10$mech<-relevel(alldata_df_2_10$mech, "path_e")
alldata_df_2_10$mech<-relevel(alldata_df_2_10$mech, "path")
alldata_df_2_10$mech<-relevel(alldata_df_2_10$mech, "dmi_e")
alldata_df_2_10$mech<-relevel(alldata_df_2_10$mech, "dmi")

### I want to do this separately for the 24 categories we have! ###
alldata_df_2_10$index<-paste(alldata_df_2_10$m, alldata_df_2_10$c, alldata_df_2_10$mech)
alldata_df_2_10$index_reps<-paste(alldata_df_2_10$m, alldata_df_2_10$c, alldata_df_2_10$mech, alldata_df_2_10$rep)
alldata_df_2_10$index<-as.factor(alldata_df_2_10$index)

alldata_df_2_10<-alldata_df_2_10[which(alldata_df_2_10$mech %in% c("dmi", "path")),]
alldata_df_2_10$rep<-as.factor(alldata_df_2_10$rep)
colours<-met.brewer(name='OKeeffe1', n=20, type='continuous') 

####Figure 4 ####
pdf(file="genomic_cline_plots_0.01_0.09_DMI_uneven.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4), mpg=c(3,3,0))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot - ERYN, rewrite the function so you can put loci number in it at some point.
for(j in 1:20){
  introgress.data<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi" & alldata_df_2_10$rep==j),c(12:14)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi" & alldata_df_2_10$rep==j),]$q
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


#### table of anovas and R2 etc.

  l1.4_genotype_fstat<-unlist(summary(aov(l1.4~rep, data_deme_2_gen_10)))[7]
  l1.4_genotype_pvalue<-unlist(summary(aov(l1.4~rep, data_deme_2_gen_10)))[9]
  R2_l1.4<-unlist(summary(aov(l1.4~as.factor(rep), data_deme_2_gen_10)))[3]/(unlist(summary(aov(l1.4~as.factor(rep), data_deme_2_gen_10)))[3]+unlist(summary(aov(l1.4~as.factor(rep), data_deme_2_gen_10)))[4])
  
  l1.10_genotype_fstat<-unlist(summary(aov(l1.10~rep, data_deme_2_gen_10)))[7]
  l1.10_genotype_pvalue<-unlist(summary(aov(l1.10~rep, data_deme_2_gen_10)))[9]
  R2_l1.10<-unlist(summary(aov(l1.10~as.factor(rep), data_deme_2_gen_10)))[3]/(unlist(summary(aov(l1.10~as.factor(rep), data_deme_2_gen_10)))[3]+unlist(summary(aov(l1.10~as.factor(rep), data_deme_2_gen_10)))[4])

  l3.4_genotype_fstat<-unlist(summary(aov(l3.4~rep, data_deme_2_gen_10)))[7]
  l3.4_genotype_pvalue<-unlist(summary(aov(l3.4~rep, data_deme_2_gen_10)))[9]
  R2_l3.4<-unlist(summary(aov(l3.4~as.factor(rep), data_deme_2_gen_10)))[3]/(unlist(summary(aov(l3.4~as.factor(rep), data_deme_2_gen_10)))[3]+unlist(summary(aov(l3.4~as.factor(rep), data_deme_2_gen_10)))[4])
  
#### just to get the triangle plots
  triangle_cols <- met.brewer("OKeeffe1", 20)
  
  triangle <- as.matrix(data_deme_2_gen_10[,c(1:8)])
  dim(triangle)
  head(triangle)
  
 # gen10 <- subset(triangle, triangle[,2]==10)
  #dim(gen10)
  #head(gen10)
  
  #gen10_deme6 <- subset(gen10, gen10[,3]==6)
  #dim(gen10_deme6)
  #head(gen10_deme6)
  
  pdf(file="triangle_0.01_0.09_DMI_uneven.pdf")
  #quartz(height=12, width=15)
  #par(mar=c(5,5,1,1), mfrow=c(4,5))
  for (i in 1:20)
  {
    rep_sub <- subset(triangle, triangle[,1]==i)
    plot(0, type="n", xlab="q", ylab="Q", xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
    segments(0, 0, 0.5, 1, lwd=3)
    segments(1, 0, 0.5, 1, lwd=3)
    for (j in 1:75) ###fewer individuals than in the regular simulations (which all have 150)
    {
      points(rep_sub[j,5], rep_sub[j,6], pch=21, bg=adjustcolor(triangle_cols[i], alpha.f=0.75), cex=2)
    }
    box(lwd=2)
  }
dev.off()
