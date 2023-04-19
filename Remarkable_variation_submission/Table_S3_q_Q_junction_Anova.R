####Script for Table S3 ####

require('stringr')

library(stringr)
library(MASS)

source("SNP_inputs.R")

nrow<-10*length(datafiles)

Fstats<-vector(length = length(unique(data_deme_6$index)))
pvalues<-list(length(unique(data_deme_6$index)))
plots<-list()

data_deme_6$index<-paste(data_deme_6$m, data_deme_6$c, data_deme_6$mech, data_deme_6$gen)

Fstats_njunct<-numeric(length = length(unique(data_deme_6$index)))
pvalues_njunct<-numeric(length(unique(data_deme_6$index)))
R2_njunction<-numeric(length(unique(data_deme_6$index)))

for(i in 1:length(unique(data_deme_6$index))){
  Fstats_njunct[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$njunct~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[7]
  pvalues_njunct[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$njunct~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[9]
  R2_njunction[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$njunct~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[3]/(unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$njunct~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[3]+unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$njunct~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[4])
}

njunct_table<-data.frame(model=as.character(unique(data_deme_6$index)),
                         Fstats=Fstats_njunct, pvalues=pvalues_njunct, Rsq=R2_njunction)

#####q score####
Fstats_q<-numeric(length = length(unique(data_deme_6$index)))
pvalues_q<-numeric(length(unique(data_deme_6$index)))
R2_q<-numeric(length(unique(data_deme_6$index)))

for(i in 1:length(unique(data_deme_6$index))){
  Fstats_q[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$q~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[7]
  pvalues_q[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$q~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[9]
  R2_q[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$q~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[3]/(unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$q~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[3]+unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$q~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[4])
}

q_table<-data.frame(model=as.character(unique(data_deme_6$index)),
                    Fstats=Fstats_q, pvalues=pvalues_q, Rsq=R2_q)


#### het (big Q)####

Fstats_het<-numeric(length = length(unique(data_deme_6$index)))
pvalues_het<-numeric(length(unique(data_deme_6$index)))
R2_het<-numeric(length(unique(data_deme_6$index)))

for(i in 1:length(unique(data_deme_6$index))){
  Fstats_het[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$het~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[7]
  pvalues_het[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$het~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[9]
  R2_het[[i]]<-unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$het~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[3]/(unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$het~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[3]+unlist(summary(aov(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$het~as.factor(data_deme_6[which(data_deme_6$index==unique(data_deme_6$index)[i]),]$rep))))[4])
}

het_table<-data.frame(model=as.character(unique(data_deme_6$index)), Fstats=Fstats_het, pvalues=pvalues_het,  Rsq=R2_het)

#save.image("/gscratch/emcfarl2/predicting_hybrids/Fstats_individuals.RData")

design.table<-data.frame(model=het_table$model,
                         migration=numeric(nrow(het_table)),
                         selection=numeric(nrow(het_table)),
                         generation=numeric(nrow(het_table)))

tmp<-matrix(as.numeric(unlist(strsplit(design.table$model, split=" "))),
            ncol=4, byrow=TRUE)
design.table$migration <- tmp[,1]
design.table$selection <- tmp[,2]
design.table$generation <- tmp[,3]

library(RColorBrewer)

mycolors1<-c(brewer.pal(9, "Purples"), "black")
mycolors2<-c(brewer.pal(9, "Oranges"), "red4")

par(mfrow=c(1,3), pty='s')
plot(log10(het_table$Fstats), log10(q_table$Fstats), type='n')
points(log10(het_table$Fstats[design.table$migration == 0.2]),
       log10(q_table$Fstats[design.table$migration == 0.2]), pch=19,
       col=mycolors1)
points(log10(het_table$Fstats[design.table$migration == 0.01]),
       log10(q_table$Fstats[design.table$migration == 0.01]), pch=19,
       col=mycolors2)

plot(log10(het_table$Fstats), log10(njunct_table$Fstats), type='n')
points(log10(het_table$Fstats[design.table$migration == 0.2]),
       log10(njunct_table$Fstats[design.table$migration == 0.2]), pch=19,
       col=mycolors1)
points(log10(het_table$Fstats[design.table$migration == 0.01]),
       log10(njunct_table$Fstats[design.table$migration == 0.01]), pch=19,
       col=mycolors2)

plot(log10(njunct_table$Fstats), log10(q_table$Fstats), type='n')
points(log10(njunct_table$Fstats[design.table$migration == 0.2]),
       log10(q_table$Fstats[design.table$migration == 0.2]), pch=19,
       col=mycolors1)
points(log10(njunct_table$Fstats[design.table$migration == 0.01]),
       log10(q_table$Fstats[design.table$migration == 0.01]), pch=19,
       col=mycolors2)

dev.print(pdf, "bivariate_ANOVA_F.pdf", width=8, height=3)

