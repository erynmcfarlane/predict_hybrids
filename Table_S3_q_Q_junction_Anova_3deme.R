####Script for Table S? ####

data_deme_2_gen_10$index<-paste(data_deme_2_gen_10$m, data_deme_2_gen_10$c, data_deme_2_gen_10$mech)

Fstats_njunct<-numeric(length = length(unique(data_deme_2_gen_10$index)))
pvalues_njunct<-numeric(length(unique(data_deme_2_gen_10$index)))

for(i in 1:length(unique(data_deme_2_gen_10$index))){
  Fstats_njunct[[i]]<-unlist(summary(aov(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$njunct~as.factor(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$rep))))[7]
  pvalues_njunct[[i]]<-unlist(summary(aov(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$njunct~as.factor(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$rep))))[9]
}

njunct_table<-data.frame(model=as.character(unique(data_deme_2_gen_10$index)),
                         Fstats=Fstats_njunct, pvalues=pvalues_njunct)

#####q score####
Fstats_q<-numeric(length = length(unique(data_deme_2_gen_10$index)))
pvalues_q<-numeric(length(unique(data_deme_2_gen_10$index)))

for(i in 1:length(unique(data_deme_2_gen_10$index))){
  Fstats_q[[i]]<-unlist(summary(aov(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$q~as.factor(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$rep))))[7]
  pvalues_q[[i]]<-unlist(summary(aov(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$q~as.factor(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$rep))))[9]
}

q_table<-data.frame(model=as.character(unique(data_deme_2_gen_10$index)),
                    Fstats=Fstats_q, pvalues=pvalues_q)


#### het (big Q)####

Fstats_het<-numeric(length = length(unique(data_deme_2_gen_10$index)))
pvalues_het<-numeric(length(unique(data_deme_2_gen_10$index)))

for(i in 1:length(unique(data_deme_2_gen_10$index))){
  Fstats_het[[i]]<-unlist(summary(aov(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$het~as.factor(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$rep))))[7]
  pvalues_het[[i]]<-unlist(summary(aov(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$het~as.factor(data_deme_2_gen_10[which(data_deme_2_gen_10$index==unique(data_deme_2_gen_10$index)[i]),]$rep))))[9]
}

het_table<-data.frame(model=as.character(unique(data_deme_2_gen_10$index)), Fstats=Fstats_het, pvalues=pvalues_het)

write.csv(het_table, file="deme_2_Q.csv", row.names = FALSE)