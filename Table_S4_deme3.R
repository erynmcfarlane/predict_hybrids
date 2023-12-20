###Script for Table S? ###

l1.4_genotype_fstat<-vector(length = length(unique(data_deme_2_10_noE$index)))
l1.4_genotype_pvalue<-vector(length = length(unique(data_deme_2_10_noE$index)))

l1.10_genotype_fstat<-vector(length = length(unique(data_deme_2_10_noE$index)))
l1.10_genotype_pvalue<-vector(length = length(unique(data_deme_2_10_noE$index)))

l3.4_genotype_fstat<-vector(length = length(unique(data_deme_2_10_noE$index)))
l3.4_genotype_pvalue<-vector(length = length(unique(data_deme_2_10_noE$index)))


for(i in 1:length(unique(data_deme_2_10_noE$index))){
  l1.4_genotype_fstat[[i]]<-unlist(summary(aov(l1.4~rep, data_deme_2_10_noE[which(data_deme_2_10_noE$index==unique(data_deme_2_10_noE$index)[i]),])))[7]
  l1.4_genotype_pvalue[[i]]<-unlist(summary(aov(l1.4~rep, data_deme_2_10_noE[which(data_deme_2_10_noE$index==unique(data_deme_2_10_noE$index)[i]),])))[9]
  
  l1.10_genotype_fstat[[i]]<-unlist(summary(aov(l1.10~rep, data_deme_2_10_noE[which(data_deme_2_10_noE$index==unique(data_deme_2_10_noE$index)[i]),])))[7]
  l1.10_genotype_pvalue[[i]]<-unlist(summary(aov(l1.10~rep, data_deme_2_10_noE[which(data_deme_2_10_noE$index==unique(data_deme_2_10_noE$index)[i]),])))[9]
  
  l3.4_genotype_fstat[[i]]<-unlist(summary(aov(l3.4~rep, data_deme_2_10_noE[which(data_deme_2_10_noE$index==unique(data_deme_2_10_noE$index)[i]),])))[7]
  l3.4_genotype_pvalue[[i]]<-unlist(summary(aov(l3.4~rep, data_deme_2_10_noE[which(data_deme_2_10_noE$index==unique(data_deme_2_10_noE$index)[i]),])))[9]
} 

Genotype.Anova<-cbind(unique(as.character(data_deme_2_10_noE$index)), l1.4_genotype_fstat, l1.4_genotype_pvalue, l1.10_genotype_fstat, l1.10_genotype_pvalue, l3.4_genotype_fstat, l3.4_genotype_pvalue)  
write.table(Genotype.Anova, file="Genotype.Anova.repfactor.Deme3.csv", col.names=T, row.names=F, quote=F, sep=',') 
