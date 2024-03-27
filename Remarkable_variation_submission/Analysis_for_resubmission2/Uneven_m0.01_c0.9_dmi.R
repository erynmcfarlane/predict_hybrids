#### between here and the next section is what's needed to run the asymetrical simulations
module load arcc/1.0
module load gcc/12.2.0
module load gsl/2.7.1
module load r

/gscratch/buerkle/data/incompatible/dfuse_src/dfuse

/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o dmi_m0.01_c0_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o dmi_m0.01_c0.2_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o dmi_m0.01_c0.9_uneven

/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o dmi_m0.2_c0_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o dmi_m0.2_c0.2_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o dmi_m0.2_c0.9_uneven

/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o path_m0.01_c0_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o path_m0.01_c0.2_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o path_m0.01_c0.9_uneven

/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o path_m0.2_c0_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o path_m0.2_c0.2_uneven
/gscratch/buerkle/data/incompatible/dfuse_src/dfuse -d 3_deme_uneven.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o path_m0.2_c0.9_uneven


#### once all of the simulations are done, move into R to do all this ####
library(tidyverse)
library(MASS)
library(data.table)
library(MetBrewer)
library(introgress)

datafiles<-list.files("/gscratch/emcfarl2/predicting_hybrids/",  pattern="*dmi_m0.01_c0.9.main", recursive=TRUE, include.dirs=TRUE)
datafiles2<-list.files("/gscratch/emcfarl2/predicting_hybrids/",  pattern="*dmi_m0.2_c0.9.main", recursive=TRUE, include.dirs=TRUE)
datafiles<-datafiles[c(1,2,3,5,6)]
datafiles2<-datafiles2[c(1,2,3,5,6)]
datafiles<-c(datafiles, datafiles2)

###this would work if I didn't have some other random stuff in there
basenames<-basename(datafiles)
m<-str_extract(basenames, "(\\d+\\.*\\d*)")
c<-str_match(basenames,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(basenames, "^([^_]+_){1}([^_])") 
u<-c("1:1","50:1","5:1","10:1","2:1", "1:1","50:1","5:1","10:1","2:1")

alldata<-list()
for(i in 1:length(datafiles)){
  alldata[[i]]<-fread(datafiles[i], sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
  alldata[[i]]$uneven<-as.factor(rep(u[i], nrow(alldata[[i]])))
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

data_deme_2_gen_10[,c(1:8, 519:522, 12, 18, 114)]->alldata_df_2_10
alldata_df_2_10$mech<-as.factor(ifelse(alldata_df_2_10$mech=="dmi_m", "dmi", alldata_df_2_10$mech))

### I want to do this separately for the 24 categories we have! ###
alldata_df_2_10$index<-paste(alldata_df_2_10$m, alldata_df_2_10$c, alldata_df_2_10$mech, alldata_df_2_10$uneven)
alldata_df_2_10$index_reps<-paste(alldata_df_2_10$m, alldata_df_2_10$c, alldata_df_2_10$mech, alldata_df_2_10$uneven, alldata_df_2_10$rep)
alldata_df_2_10$index<-as.factor(alldata_df_2_10$index)

alldata_df_2_10<-alldata_df_2_10[which(alldata_df_2_10$mech %in% c("dmi", "path")),]
alldata_df_2_10$rep<-as.factor(alldata_df_2_10$rep)
colours<-met.brewer(name='OKeeffe1', n=20, type='continuous') 

####Figure 4 #### 1:1
pdf(file="genomic_cline_plots_0.01_0.09_DMI_uneven1:1.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot)
for(j in 1:20){
  introgress.data<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1" & alldata_df_2_10$rep==j),c(13:15)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1" & alldata_df_2_10$rep==j),]$q
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

####Figure 4 #### 2:1
pdf(file="genomic_cline_plots_0.01_0.09_DMI_uneven2:1.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot)
for(j in 1:20){
  introgress.data<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1" & alldata_df_2_10$rep==j),c(13:15)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1" & alldata_df_2_10$rep==j),]$q
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

####Figure 4 #### 5:1
pdf(file="genomic_cline_plots_0.01_0.09_DMI_uneven5:1.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot)
for(j in 1:20){
  introgress.data<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1" & alldata_df_2_10$rep==j),c(13:15)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1" & alldata_df_2_10$rep==j),]$q
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

####Figure 4 #### 10:1
pdf(file="genomic_cline_plots_0.01_0.09_DMI_uneven10:1.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot)
for(j in 1:20){
  introgress.data<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1" & alldata_df_2_10$rep==j),c(13:15)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1" & alldata_df_2_10$rep==j),]$q
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

####Figure 4 #### 50:1
pdf(file="genomic_cline_plots_0.01_0.09_DMI_uneven50:1.pdf", width=30, height=10)

par(mfrow=c(1,3), mar=c(5,5,1,1), oma=c(5,5,4,4))
genomic.clines.reps<-list()
Fitted.AA<-list()
Fitted.Aa<-list()
### do this three times, and change the loci in the genomic cline plot)
for(j in 1:20){
  introgress.data<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1" & alldata_df_2_10$rep==j),c(13:15)]
  chromosome<-(as.numeric(unlist(str_extract(colnames(introgress.data),"[[:digit:]]+\\."))))
  hi.index<-alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1" & alldata_df_2_10$rep==j),]$q
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

  l1.4_genotype_fstat_0.01_50<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[7]
    R2_l1.4_0.01_50<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[4])
  
  l1.10_genotype_fstat_0.01_50<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[7]
   R2_l1.10_0.01_50<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[4])
  
  l3.4_genotype_fstat_0.01_50<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[7]
  R2_l3.4_0.01_50<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 50:1"),])))[4])
  
  l1.4_genotype_fstat_0.2_50<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[7]
   R2_l1.4_0.2_50<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[4])
  
  l1.10_genotype_fstat_0.2_50<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[7]
  R2_l1.10_0.2_50<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[4])
  
  l3.4_genotype_fstat_0.2_50<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[7]
  R2_l3.4_0.2_50<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 50:1"),])))[4])
  
  
  l1.4_genotype_fstat_0.01_10<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[7]
  R2_l1.4_0.01_10<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[4])
  
  l1.10_genotype_fstat_0.01_10<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[7]
  R2_l1.10_0.01_10<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[4])
  
  l3.4_genotype_fstat_0.01_10<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[7]
  R2_l3.4_0.01_10<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 10:1"),])))[4])
  
  l1.4_genotype_fstat_0.2_10<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[7]
  R2_l1.4_0.2_10<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[4])
  
  l1.10_genotype_fstat_0.2_10<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[7]
  R2_l1.10_0.2_10<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[4])
  
  l3.4_genotype_fstat_0.2_10<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[7]
  R2_l3.4_0.2_10<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 10:1"),])))[4])
  
  
  l1.4_genotype_fstat_0.01_5<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[7]
  R2_l1.4_0.01_5<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[4])
  
  l1.10_genotype_fstat_0.01_5<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[7]
  R2_l1.10_0.01_5<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[4])
  
  l3.4_genotype_fstat_0.01_5<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[7]
  R2_l3.4_0.01_5<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 5:1"),])))[4])
  
  l1.4_genotype_fstat_0.2_5<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[7]
  R2_l1.4_0.2_5<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[4])
  
  l1.10_genotype_fstat_0.2_5<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[7]
  R2_l1.10_0.2_5<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[4])
  
  l3.4_genotype_fstat_0.2_5<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[7]
  R2_l3.4_0.2_5<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 5:1"),])))[4])
  
  
  l1.4_genotype_fstat_0.01_2<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[7]
  R2_l1.4_0.01_2<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[4])
  
  l1.10_genotype_fstat_0.01_2<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[7]
  R2_l1.10_0.01_2<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[4])
  
  l3.4_genotype_fstat_0.01_2<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[7]
  R2_l3.4_0.01_2<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 2:1"),])))[4])
  
  l1.4_genotype_fstat_0.2_2<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[7]
  R2_l1.4_0.2_2<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[4])
  
  l1.10_genotype_fstat_0.2_2<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[7]
  R2_l1.10_0.2_2<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[4])
  
  l3.4_genotype_fstat_0.2_2<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[7]
  R2_l3.4_0.2_2<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 2:1"),])))[4])
  
  
  l1.4_genotype_fstat_0.01_1<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[7]
  R2_l1.4_0.01_1<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[4])
  
  l1.10_genotype_fstat_0.01_1<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[7]
  R2_l1.10_0.01_1<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[4])
  
  l3.4_genotype_fstat_0.01_1<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[7]
  R2_l3.4_0.01_1<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.01 0.9 dmi 1:1"),])))[4])
  
  l1.4_genotype_fstat_0.2_1<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[7]
  R2_l1.4_0.2_1<-unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[4])
  
  l1.10_genotype_fstat_0.2_1<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[7]
  R2_l1.10_0.2_1<-unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[4])
  
  l3.4_genotype_fstat_0.2_1<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[7]
  R2_l3.4_0.2_1<-unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_2_10[which(alldata_df_2_10$index=="0.2 0.9 dmi 1:1"),])))[4])
  
SNP<-c(rep("1.4", 10), rep("1.10", 10), rep("3.4", 10))
m_reps<-c(rep(c(rep(0.01, 5), rep(0.2, 5)), 3))
asymmetry<-rep(c(50, 10, 5, 2, 1), 6)
plot_labels<-rep(c(0.5, 1.5, 2.5, 3.5, 4.5), 6)
fstat<-c(l1.4_genotype_fstat_0.01_50, l1.4_genotype_fstat_0.01_10, l1.4_genotype_fstat_0.01_5, l1.4_genotype_fstat_0.01_2, l1.4_genotype_fstat_0.01_1, l1.4_genotype_fstat_0.2_50, l1.4_genotype_fstat_0.2_10, l1.4_genotype_fstat_0.2_5, l1.4_genotype_fstat_0.2_2, l1.4_genotype_fstat_0.2_1, l1.10_genotype_fstat_0.01_50, l1.10_genotype_fstat_0.01_10, l1.10_genotype_fstat_0.01_5, l1.10_genotype_fstat_0.01_2, l1.10_genotype_fstat_0.01_1, l1.10_genotype_fstat_0.2_50, l1.10_genotype_fstat_0.2_10, l1.10_genotype_fstat_0.2_5, l1.10_genotype_fstat_0.2_2, l1.10_genotype_fstat_0.2_1, l3.4_genotype_fstat_0.01_50, l3.4_genotype_fstat_0.01_10, l3.4_genotype_fstat_0.01_5, l3.4_genotype_fstat_0.01_2, l3.4_genotype_fstat_0.01_1, l3.4_genotype_fstat_0.2_50, l3.4_genotype_fstat_0.2_10, l3.4_genotype_fstat_0.2_5, l3.4_genotype_fstat_0.2_2, l3.4_genotype_fstat_0.2_1)
R2<-c(R2_l1.4_0.01_50,  R2_l1.4_0.01_10,  R2_l1.4_0.01_5,  R2_l1.4_0.01_2,  R2_l1.4_0.01_1, R2_l1.4_0.2_50,  R2_l1.4_0.2_10,  R2_l1.4_0.2_5,  R2_l1.4_0.2_2,  R2_l1.4_0.2_1, R2_l1.10_0.01_50,  R2_l1.10_0.01_10,  R2_l1.10_0.01_5,  R2_l1.10_0.01_2,  R2_l1.10_0.01_1,R2_l1.10_0.2_50,  R2_l1.10_0.2_10,  R2_l1.10_0.2_5,  R2_l1.10_0.2_2,  R2_l1.10_0.2_1, R2_l3.4_0.01_50,  R2_l3.4_0.01_10,  R2_l3.4_0.01_5,  R2_l3.4_0.01_2,  R2_l3.4_0.01_1,R2_l3.4_0.2_50,  R2_l3.4_0.2_10,  R2_l3.4_0.2_5,  R2_l3.4_0.2_2,  R2_l3.4_0.2_1)

cbind.data.frame(SNP, m_reps, asymmetry, fstat, R2, plot_labels)->asym_compare

pdf(file="locus_specific_Fstat_R2_asym_compare.pdf", width=6, height=9)

par(mfcol=c(3,2), mar=c(5,5,1,1))
plot(0, type="n", xlab="degree of asymmetry", ylab="Fstat - SNP 1.4", xlim=c(0, 5), ylim=c(0, 60), las=1, cex.lab=1.5, cex.axis=1.25, xaxt="n")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5), label=c("1:50", "1:10", "1:5", "1:2", "1:1"), cex.axis=1.15)
lines(fstat~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.4" & asym_compare$m_reps==0.01),], pch=21,bg="red", cex=1.5)
lines(fstat~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.4" & asym_compare$m_reps==0.2),], pch=21,bg="blue", cex=1.5)

plot(0, type="n", xlab="degree of asymmetry", ylab="Fstat - SNP 1.10", xlim=c(0, 5), ylim=c(0, 60), las=1, cex.lab=1.5, cex.axis=1.25, xaxt="n")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5), label=c("1:50", "1:10", "1:5", "1:2", "1:1"), cex.axis=1.15)
lines(fstat~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.10" & asym_compare$m_reps==0.01),], pch=21,bg="red", cex=1.5)
lines(fstat~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.10" & asym_compare$m_reps==0.2),], pch=21,bg="blue", cex=1.5)

plot(0, type="n", xlab="degree of asymmetry", ylab="Fstat - SNP 3.4", xlim=c(0, 5), ylim=c(0, 60), las=1, cex.lab=1.5, cex.axis=1.25, xaxt="n")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5), label=c("1:50", "1:10", "1:5", "1:2", "1:1"), cex.axis=1.15)
lines(fstat~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="3.4" & asym_compare$m_reps==0.01),], pch=21,bg="red", cex=1.5)
lines(fstat~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="3.4" & asym_compare$m_reps==0.2),], pch=21,bg="blue", cex=1.5)


plot(0, type="n", xlab="degree of asymmetry", ylab="R2", xlim=c(0, 5), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25, xaxt="n")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5), label=c("1:50", "1:10", "1:5", "1:2", "1:1"), cex.axis=1.15)
lines(R2~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.4" & asym_compare$m_reps==0.01),], pch=21,bg="red", cex=1.5)
lines(R2~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.4" & asym_compare$m_reps==0.2),], pch=21,bg="blue", cex=1.5)

plot(0, type="n", xlab="degree of asymmetry", ylab="R2", xlim=c(0, 5), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25, xaxt="n")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5), label=c("1:50", "1:10", "1:5", "1:2", "1:1"), cex.axis=1.15)
lines(R2~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.10" & asym_compare$m_reps==0.01),], pch=21,bg="red", cex=1.5)
lines(R2~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="1.10" & asym_compare$m_reps==0.2),], pch=21,bg="blue", cex=1.5)

plot(0, type="n", xlab="degree of asymmetry", ylab="R2", xlim=c(0, 5), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25, xaxt="n")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5), label=c("1:50", "1:10", "1:5", "1:2", "1:1"), cex.axis=1.15)
lines(R2~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="3.4" & asym_compare$m_reps==0.01),], pch=21,bg="red", cex=1.5)
lines(R2~plot_labels, type='o', data=asym_compare[which(asym_compare$SNP=="3.4" & asym_compare$m_reps==0.2),], pch=21,bg="blue", cex=1.5)

dev.off()

