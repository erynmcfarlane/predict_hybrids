### depends on SNP_inputs.R

alldata_df[,c(1:8,519:521)]->alldata_df
alldata_df[which(alldata_df$gen==10),]->alldata_df_gen_10
alldata_df_gen_10[which(alldata_df_gen_10$mech=="path_m" |alldata_df_gen_10$mech=="dmi_m"),]->alldata_df_gen_10_noe


alldata_df_gen_10_noe$index<-paste(alldata_df_gen_10_noe$m, alldata_df_gen_10_noe$c, alldata_df_gen_10_noe$mech, alldata_df_gen_10_noe$rep)

matrix<-matrix(data=NA, nrow=length(unique(alldata_df_gen_10_noe$index)), ncol=3)

for(i in 1:length(unique(alldata_df_gen_10_noe$index))){
  index<-data.frame(alldata_df_gen_10_noe[which(alldata_df_gen_10_noe$index==unique(alldata_df_gen_10_noe$index)[i],)])
  deme5<-ifelse(index[which(index$deme==5),]$q==1 | index[which(index$deme==5),]$q==0, "parental", "hybrid")
  deme6<-ifelse(index[which(index$deme==6),]$q==1 | index[which(index$deme==6),]$q==0, "parental", "hybrid")
  deme7<-ifelse(index[which(index$deme==7),]$q==1 | index[which(index$deme==7),]$q==0, "parental", "hybrid")
  matrix[i,1]<-sum(deme5=="hybrid")
  matrix[i,2]<-sum(deme6=="hybrid")
  matrix[i,3]<-sum(deme7=="hybrid")
  print(unique(alldata_df_gen_10_noe$index)[i])
}


sum(matrix[,2]>=matrix[,1])
sum(matrix[,2]>=matrix[,3])


#### 100 generations ###

alldata_df[which(alldata_df$gen==100),]->alldata_df_gen_100
alldata_df_gen_100[which(alldata_df_gen_100$mech=="path_m" |alldata_df_gen_100$mech=="dmi_m"),]->alldata_df_gen_100_noe


alldata_df_gen_100_noe$index<-paste(alldata_df_gen_100_noe$m, alldata_df_gen_100_noe$c, alldata_df_gen_100_noe$mech, alldata_df_gen_100_noe$rep)

matrix_100<-matrix(data=NA, nrow=length(unique(alldata_df_gen_100_noe$index)), ncol=3)

for(i in 1:length(unique(alldata_df_gen_100_noe$index))){
  index<-data.frame(alldata_df_gen_100_noe[which(alldata_df_gen_100_noe$index==unique(alldata_df_gen_100_noe$index)[i],)])
  deme5<-ifelse(index[which(index$deme==5),]$q==1 | index[which(index$deme==5),]$q==0, "parental", "hybrid")
  deme6<-ifelse(index[which(index$deme==6),]$q==1 | index[which(index$deme==6),]$q==0, "parental", "hybrid")
  deme7<-ifelse(index[which(index$deme==7),]$q==1 | index[which(index$deme==7),]$q==0, "parental", "hybrid")
  matrix_100[i,1]<-sum(deme5=="hybrid")
  matrix_100[i,2]<-sum(deme6=="hybrid")
  matrix_100[i,3]<-sum(deme7=="hybrid")
  print(unique(alldata_df_gen_100_noe$index)[i])
}


sum(matrix_100[,2]>=matrix_100[,1])
sum(matrix_100[,2]>=matrix_100[,3])