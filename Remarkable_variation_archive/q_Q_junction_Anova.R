require('stringr')

library(stringr)
library(MASS)
#library(fitdistrplus)
##if using local computer
# datafiles<-list.files("~/Google Drive/Replicate Hybrid zone review/Predicting_Hybrids_analysis/hybrid_sims/deme11_all" , pattern="*first8cols.txt.gz")
##if using teton
datafiles<-list.files("/gscratch/buerkle/data/incompatible/runs/deme11",  pattern="*main")
m<-str_extract(datafiles, "(\\d+\\.*\\d*)")
c<-str_match(datafiles,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
mech<-str_extract(datafiles, "^([^_]+_){1}([^_])") ### where e means there's an environmental interaction and m means there isn't

alldata<-list()

# setwd("/gscratch/buerkle/data/incompatible/runs/deme11")

for(i in 1:length(datafiles)){
  alldata[[i]]<-read.table(gzfile(datafiles[i]), sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}

setwd("/gscratch/emcfar2/predicting_hybrids")
alldata_df<-do.call(rbind.data.frame, alldata)
summary(alldata_df)
### I'm not at all sure that I need this, so maybe take it out if it's not needed. 
alldata_df$gen<-as.factor(alldata_df$gen)
#need to take the mean # of junctions for 19 replicates, for each generation, for each deme, for each m, for each c, and ask how the individual number of junctions for the other replicate fits the distribution
#do this for each replicate, through all 20
nrow<-10*length(datafiles)

Fstats<-vector(length = length(unique(alldata_df_6$index)))
pvalues<-list(length(unique(alldata_df_6$index)))
plots<-list()
#only going to use middle, hybridizing demes 5 and 6
#two reasons for this - 1), this is really what we're most interested in, and 2) we were getting -INFs when the mean njuncts was 0. 
alldata_df_6<-subset(alldata_df, deme %in% c(6))
alldata_df_6$index<-paste(alldata_df_6$m, alldata_df_6$c, alldata_df_6$gen, alldata_df_6$deme)

Fstats_njunct<-numeric(length = length(unique(alldata_df_6$index)))
pvalues_njunct<-numeric(length(unique(alldata_df_6$index)))

for(i in 1:length(unique(alldata_df_6$index))){
  Fstats_njunct[[i]]<-unlist(summary(aov(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$njunct~as.factor(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$rep))))[7]
  pvalues_njunct[[i]]<-unlist(summary(aov(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$njunct~as.factor(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$rep))))[9]
}

njunct_table<-data.frame(model=as.character(unique(alldata_df_6$index)),
                         Fstats=Fstats_njunct, pvalues=pvalues_njunct)

#####q score####
Fstats_q<-numeric(length = length(unique(alldata_df_6$index)))
pvalues_q<-numeric(length(unique(alldata_df_6$index)))

for(i in 1:length(unique(alldata_df_6$index))){
  Fstats_q[[i]]<-unlist(summary(aov(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$q~as.factor(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$rep))))[7]
  pvalues_q[[i]]<-unlist(summary(aov(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$q~as.factor(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$rep))))[9]
}

q_table<-data.frame(model=as.character(unique(alldata_df_6$index)),
                    Fstats=Fstats_q, pvalues=pvalues_q)


#### het (big Q)####

Fstats_het<-numeric(length = length(unique(alldata_df_6$index)))
pvalues_het<-numeric(length(unique(alldata_df_6$index)))

for(i in 1:length(unique(alldata_df_6$index))){
  Fstats_het[[i]]<-unlist(summary(aov(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$het~as.factor(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$rep))))[7]
  pvalues_het[[i]]<-unlist(summary(aov(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$het~as.factor(alldata_df_6[which(alldata_df_6$index==unique(alldata_df_6$index)[i]),]$rep))))[9]
}

het_table<-data.frame(model=as.character(unique(alldata_df_6$index)), Fstats=Fstats_het, pvalues=pvalues_het)

save.image("/gscratch/emcfarl2/predicting_hybrids/Fstats_individuals.RData")

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

