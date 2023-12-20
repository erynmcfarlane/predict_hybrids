### These are all the minor ancestry analyses, giving us Table S5

###Must run SNP_inputs.R - then this script will work.
source("SNP_inputs.R")

data_deme_6_gen_10[,c(1:8, 519:521, 12, 18, 114)]->alldata_df_6_10 ### This gives the 3 loci of interest that we use for TabS4, TabS5, Fig4, FigS8-FigS10
data_long<-gather(alldata_df_6_10, snp, genotype, l1.4:l3.4, factor_key=TRUE)

source('summarySE.R')

data_long$mech<-relevel(data_long$mech, "path_e")
data_long$mech<-relevel(data_long$mech, "path_m")
data_long$mech<-relevel(data_long$mech, "dmi_e")
data_long$mech<-relevel(data_long$mech, "dmi_m")
data_long$mech<-ifelse(data_long$mech=="dmi_m", "dmi", ifelse(data_long$mech=="path_m", "path",ifelse(data_long$mech=="path_e", 'path_e', "dmi_e")))

data_long$snp_num<-(as.numeric(unlist(str_extract(as.factor(data_long$snp),"[[:digit:]]+\\.*[[:digit:]]*"))))
summary(data_long$snp_num)
data_long$q<-data_long$genotype/2
data_long$Q<-ifelse(data_long$genotype==1, 1, 0) 
data_long_noE<-data_long[which(data_long$mech %in% c("dmi", "path")),]
summaries_q<-summarySE(data_long_noE, measurevar='q', groupvars=c("m", "c", "mech", 'rep'), na.rm=FALSE, conf.interval=.95)
summaries_q$partial_index<-paste(summaries_q$m, summaries_q$c, summaries_q$mech, summaries_q$rep)

data_long_noE$partial_index<-paste(data_long_noE$m, data_long_noE$c, data_long_noE$mech, data_long_noE$rep)

merge(data_long_noE, summaries_q[,c(10, 6)], by='partial_index')->data_long_merged
data_long_merged$rep<-as.factor(data_long_merged$rep)
data_long_merged$minorancestry<-ifelse(data_long_merged$q.y<=0.5, data_long_merged$q.x, 1-data_long_merged$q.x)
colours<-met.brewer(name='OKeeffe1', n=20, type='continuous') 

### want to do 1) correlations between populations and 2) anovas between populations for each SNP
### the schumer papers do this based on haplotypes, I think, not genotypes. I'm going to set it up for genotypes, but then I can do windows? Except we don't have a linkage map...
l1.4_minorancestry_fstat<-vector(length = length(unique(data_long_merged$scenario)))
l1.4_minorancestry_pvalue<-vector(length = length(unique(data_long_merged$scenario)))
correlation_1.4<-vector(length = length(unique(data_long_merged$scenario)))

l1.10_minorancestry_fstat<-vector(length = length(unique(data_long_merged$scenario)))
l1.10_minorancestry_pvalue<-vector(length = length(unique(data_long_merged$scenario)))
correlation_1.10<-vector(length = length(unique(data_long_merged$scenario)))

l3.4_minorancestry_fstat<-vector(length = length(unique(data_long_merged$scenario)))
l3.4_minorancestry_pvalue<-vector(length = length(unique(data_long_merged$scenario)))
correlation_3.4<-vector(length = length(unique(data_long_merged$scenario)))

for(i in 1:length(unique(data_long_merged$scenario))){
  l1.4_minorancestry_fstat[[i]]<-unlist(summary(aov(minorancestry~rep, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4),])))[7]
  l1.4_minorancestry_pvalue[[i]]<-unlist(summary(aov(minorancestry~rep,data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4),])))[9]
  
  l1.10_minorancestry_fstat[[i]]<-unlist(summary(aov(minorancestry~rep, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10),])))[7]
  l1.10_minorancestry_pvalue[[i]]<-unlist(summary(aov(minorancestry~rep,data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10),])))[9]
  
  l3.4_minorancestry_fstat[[i]]<-unlist(summary(aov(minorancestry~rep, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4),])))[7]
  l3.4_minorancestry_pvalue[[i]]<-unlist(summary(aov(minorancestry~rep,data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4),])))[9]
  
  ###let's also do a correlation between the ancestries for rep 1 and 2 for each scenario for each snp
  correlation_1.4[[i]]<-cor.test(data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4 & data_long_merged$rep==1),]$minorancestry, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.4 & data_long_merged$rep==2),]$minorancestry)$estimate
  correlation_1.10[[i]]<-cor.test(data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10 & data_long_merged$rep==1),]$minorancestry, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==1.10 & data_long_merged$rep==2),]$minorancestry)$estimate
  correlation_3.4[[i]]<-cor.test(data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4 & data_long_merged$rep==1),]$minorancestry, data_long_merged[which(data_long_merged$scenario==unique(data_long_merged$scenario)[i] & data_long_merged$snp_num==3.4 & data_long_merged$rep==2),]$minorancestry)$estimate
  

  }

minorancestry.Anova<-cbind(unique(as.character(data_long_merged$scenario)), l1.4_minorancestry_fstat, l1.4_minorancestry_pvalue, correlation_1.4, l1.10_minorancestry_fstat, l1.10_minorancestry_pvalue, correlation_1.10, l3.4_minorancestry_fstat, l3.4_minorancestry_pvalue, correlation_3.4)  
write.table(minorancestry.Anova, file="minorancestry.Anova.csv", col.names=T, row.names=F, quote=F, sep=',') 
