#### Script that gives Figure 3, different inputs gives Figures S3 A and B

###Must run SNP_inputs.R - then this script will work.
source("SNP_inputs.R")


data_long<-gather(data_deme_6_gen_10, snp, genotype, l1.1:l10.51, factor_key=TRUE)

source('summarySE.R')

data_long$mech<-relevel(data_long$mech, "path_e")
data_long$mech<-relevel(data_long$mech, "path_m")
data_long$mech<-relevel(data_long$mech, "dmi_e")
data_long$mech<-relevel(data_long$mech, "dmi_m")
data_long$mech<-ifelse(data_long$mech=="dmi_m", "dmi", ifelse(data_long$mech=="path_m", "path",ifelse(data_long$mech=="path_e", 'path_e', "dmi_e")))

data_long$snp_num<-(as.numeric(unlist(str_extract(as.factor(data_long$snp),"[[:digit:]]+\\.*[[:digit:]]*"))))
summary(data_long$snp_num)
data_long$q<-data_long$genotype/2
data_long$Q<-ifelse(data_long$genotype==1, 1, 0) ### is this right? There are really only two locus-specific options, right?

data_long_noE<-data_long[which(data_long$mech %in% c("dmi", "path")),]
summaries_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'rep','snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_q$partial_index<-paste(summaries_q$m, summaries_q$c, summaries_q$mech, summaries_q$snp_num)
summaries_q$index_nosnp<-paste(summaries_q$m, summaries_q$c, summaries_q$mech)
summaries_mean_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_q$partial_index<-paste(summaries_mean_q$m, summaries_mean_q$c, summaries_mean_q$mech, summaries_mean_q$snp_num)
summaries_mean_q<-summaries_mean_q[,c(10, 6)]
names(summaries_mean_q)<-c("partial_index", "mean_q")
summaries_q<-merge(summaries_q,summaries_mean_q, by='partial_index')

summaries_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'rep','snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_Q$partial_index<-paste(summaries_Q$m, summaries_Q$c, summaries_Q$mech, summaries_Q$snp_num)
summaries_Q$index_nosnp<-paste(summaries_Q$m, summaries_Q$c, summaries_Q$mech)
summaries_mean_Q<-summarySE(data_long_noE, measurevar='Q', groupvars=c("m", "c", "mech", 'snp_num'), na.rm=FALSE, conf.interval=.95)
summaries_mean_Q$partial_index<-paste(summaries_mean_Q$m, summaries_mean_Q$c, summaries_mean_Q$mech, summaries_mean_Q$snp_num)
summaries_mean_Q<-summaries_mean_Q[,c(10, 6)]
names(summaries_mean_Q)<-c("partial_index", "mean_Q")
summaries_Q<-merge(summaries_Q,summaries_mean_Q, by='partial_index')

summaries_q<-summaries_q[which(summaries_q$snp_num<4),] ### to only use the chromosomes I'm gonna plot
summaries_Q<-summaries_Q[which(summaries_Q$snp_num<4),]

write.csv(summaries_q, file="summaries_q.csv", quote=FALSE, row.names=FALSE)
write.csv(summaries_Q, file="summaries_Q12.csv", quote=FALSE, row.names=FALSE)


### ------------------ beginning of plotting

summaries_q<-read.csv("summaries_q.csv")
summaries_Q<-read.csv("summaries_Q12.csv")

library(MetBrewer)
#library(patchwork)
colours <- met.brewer(name="OKeeffe1", n=20, type="continuous")


###admix --------------------------------------------------------
pdf(file="plotsum_admix_10.pdf", width=10, height=6)
layout(matrix(c(13, 16:21, 14, seq(1,11,2), 15, seq(2,12,2)), nrow=3, byrow=TRUE),
       widths=c(3,rep(5,6)), heights=c(1, 5, 5))
par(mar=c(4,2,0.1,0.1))
for(i in 1:length(unique(summaries_q$index_nosnp))){
  plot(0, type="l", xlab="", ylab="", ylim=c(0,1), xlim=c(1.1,3.91), axes=F)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col='light gray')
  box()
  axis(1, at=c(1.5,2.5,3.5), labels=c(1,2, 3), tick=FALSE, cex.axis=1.5)
  axis(1, at=c(1.1,2,3,3.95), labels=c("", "", "", ""))
  summaries_q_temp <- summaries_q[which(summaries_q$index_nosnp ==
                                        unique(summaries_q$index_nosnp)[i]),]
  for(j in 1:20){
    ## cab: I broke the connection between chromosome 1 and 2 markers
    ### sem: since the path models include chromosome 2, we're going to take this through all of chromosome 3, even if they're no longer consistent with Doro's figure
    lines(summaries_q_temp[which(summaries_q_temp$rep==j),6][1:46] ,
          summaries_q_temp[which(summaries_q_temp$rep==j),8][1:46], col=colours[j])
    lines(summaries_q_temp[which(summaries_q_temp$rep==j),6][47:92] ,
          summaries_q_temp[which(summaries_q_temp$rep==j),8][47:92], col=colours[j])
    lines(summaries_q_temp[which(summaries_q_temp$rep==j),6][93:138] ,
          summaries_q_temp[which(summaries_q_temp$rep==j),8][93:138], col=colours[j])
    
    lines(summaries_q_temp[which(summaries_q_temp$rep==j),6][1:46],
          summaries_q_temp[which(summaries_q_temp$rep==j),13][1:46], col='black')
    lines(summaries_q_temp[which(summaries_q_temp$rep==j),6][47:92],
          summaries_q_temp[which(summaries_q_temp$rep==j),13][47:92], col='black')
    lines(summaries_q_temp[which(summaries_q_temp$rep==j),6][93:138],
          summaries_q_temp[which(summaries_q_temp$rep==j),13][93:138], col='black')
  }
  if(i %in% 1:2){
    axis(2, at=c(0,0.5,1), labels=c(0, "", 1), cex.axis=1.5)
  }

  if (i %in% c(1,3,5,7,9,11)){
    mtext(paste0("m = ",
                 summaries_q[which(summaries_q$index==unique(summaries_q$index)[i]),]$m[i]),
          line=2, cex=1.3)
    mtext(paste0("c = ", summaries_q[which(summaries_q$index==unique(summaries_q$index)[i]),]$c[i]), line=0.25, cex=1.3)
  } else{
    mtext('Chromosome', side=1, line = 3)
  }
}
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) ## plot 13
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) ## plot 14
text(0.7, 0.5, "Admixture proportion (q)", srt=90, cex=1.5)
text(0.2, 0.5, "BDMI", cex=2, srt=90)
mtext("A)", side=3, line=1, cex=2.5)
plot(0:1, 0:1,  type="n", xlab="", ylab="", axes=FALSE) # plot 15
text(0.7, 0.5, "Admixture proportion (q)", srt=90, cex=1.5)
text(0.2, 0.5, "path", cex=2, srt=90)
## did not plot 16:21, but instead simply let mtext bleed in
dev.off()

### Intersource Ancestry ### ---------------------------------------
pdf(file="plotsum_intersource_10.pdf", width=10, height=6)
layout(matrix(c(13, 16:21, 14, seq(1,11,2), 15, seq(2,12,2)), nrow=3, byrow=TRUE),
       widths=c(3,rep(5,6)), heights=c(1, 5, 5))
par(mar=c(4,2,0.1,0.1))
for(i in 1:length(unique(summaries_Q$index_nosnp))){
  plot(0, type="l", xlab="", ylab="", ylim=c(0,1), xlim=c(1.1,3.91), axes=F)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col='light gray')
  box()
  axis(1, at=c(1.5,2.5,3.5), labels=c(1,2, 3), tick=FALSE, cex.axis=1.5)
  axis(1, at=c(1.1,2,3,3.95), labels=c("", "", "", ""))
  summaries_Q_temp<-summaries_Q[which(summaries_Q$index_nosnp == unique(summaries_Q$index_nosnp)[i]),]
  for(j in 1:20){
    ## cab: I broke the connection between chromosome 1 and 2 markers
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6][1:46] ,
          summaries_Q_temp[which(summaries_Q_temp$rep==j),8][1:46], col=colours[j])
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6][47:92] ,
          summaries_Q_temp[which(summaries_Q_temp$rep==j),8][47:92], col=colours[j])
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6][93:138] ,
          summaries_Q_temp[which(summaries_Q_temp$rep==j),8][93:138], col=colours[j])
    
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6][1:46],
          summaries_Q_temp[which(summaries_Q_temp$rep==j),13][1:46], col='black')
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6][47:92],
          summaries_Q_temp[which(summaries_Q_temp$rep==j),13][47:92], col='black')
    lines(summaries_Q_temp[which(summaries_Q_temp$rep==j),6][93:138],
          summaries_Q_temp[which(summaries_Q_temp$rep==j),13][93:138], col='black')
  }
  if(i %in% 1:2){
    axis(2, at=c(0,0.5,1), labels=c(0, "", 1), cex.axis=1.5)
  }

  if (i %in% c(1,3,5,7,9,11)){
    mtext(paste0("m = ",
                 summaries_Q[which(summaries_Q$index==unique(summaries_Q$index)[i]),]$m[i]),
          line=2, cex=1.3)
    mtext(paste0("c = ", summaries_Q[which(summaries_Q$index==unique(summaries_Q$index)[i]),]$c[i]), line=0.25, cex=1.3)
  } else{
    mtext('Chromosome', side=1, line = 3)
  }
}
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) ## plot 13
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) ## plot 14
text(0.7, 0.5, "Interspecific ancestry (Q12)", srt=90, cex=1.5)
text(0.2, 0.5, "BDMI", cex=2, srt=90)
mtext("B)", side=3, line=1, cex=2.5)
plot(0:1, 0:1,  type="n", xlab="", ylab="", axes=FALSE) # plot 15
text(0.7, 0.5, "Interspecific ancestry (Q12)", srt=90, cex=1.5)
text(0.2, 0.5, "path", cex=2, srt=90)
## did not plot 16:21, but instead simply let mtext bleed in
dev.off()
