###Script for Table S4 ###

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

l1.4_genotype_fstat<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l1.4_genotype_pvalue<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l1.4_genotype_R2<-vector(length = length(unique(alldata_df_6_10_noE$index)))

l1.10_genotype_fstat<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l1.10_genotype_pvalue<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l1.10_genotype_R2<-vector(length = length(unique(alldata_df_6_10_noE$index)))

l3.4_genotype_fstat<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l3.4_genotype_pvalue<-vector(length = length(unique(alldata_df_6_10_noE$index)))
l3.4_genotype_R2<-vector(length = length(unique(alldata_df_6_10_noE$index)))


for(i in 1:length(unique(alldata_df_6_10_noE$index))){
  l1.4_genotype_fstat[[i]]<-unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[7]
  l1.4_genotype_pvalue[[i]]<-unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[9]
  l1.4_genotype_R2[[i]]<-unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[3]/(unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[3]+unlist(summary(aov(l1.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[4])
  
  l1.10_genotype_fstat[[i]]<-unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[7]
  l1.10_genotype_pvalue[[i]]<-unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[9]
  l1.10_genotype_R2[[i]]<-unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[3]/(unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[3]+unlist(summary(aov(l1.10~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[4])

  l3.4_genotype_fstat[[i]]<-unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[7]
  l3.4_genotype_pvalue[[i]]<-unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[9]
  l3.4_genotype_R2[[i]]<-unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[3]/(unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[3]+unlist(summary(aov(l3.4~rep, alldata_df_6_10_noE[which(alldata_df_6_10_noE$index==unique(alldata_df_6_10_noE$index)[i]),])))[4])
  
} 

Genotype.Anova<-cbind(unique(as.character(alldata_df_6_10_noE$index)), l1.4_genotype_fstat, l1.4_genotype_pvalue,l1.4_genotype_R2,l1.10_genotype_fstat, l1.10_genotype_pvalue,l1.10_genotype_R2, l3.4_genotype_fstat, l3.4_genotype_pvalue, l3.4_genotype_R2)  
write.table(Genotype.Anova, file="Genotype.Anova.repfactor.csv", col.names=T, row.names=F, quote=F, sep=',') 
