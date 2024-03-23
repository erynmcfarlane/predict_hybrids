###stats on stats (but not in a bad way) for uneven parental simulations

read.csv("./Analysis_for_resubmission2/Summary_stats_long.csv")->summary_stats
head(summary_stats)
paste0(summary_stats$Architecture, summary_stats$m, summary_stats$c, summary_stats$Gen)->summary_stats$scenario


### need to switch the pvalues because lots are "<0.001"

###this doesn't work yet, need to figure out what I want from this, if I want to make a plot
for(i in 1:length(summary_stats$scenario)){
  Fstat<-lm(summary_stats[which(summary_stats$scenario==unique(summary_stats$scenario[i])),]$qF~summary_stats[which(summary_stats$scenario==unique(summary_stats$scenario[i])),]$parental_contribution), 
  pvalue<-lm(summary_stats[which(summary_stats$scenario==unique(summary_stats$scenario[i])),]$q.P~summary_stats[which(summary_stats$scenario==unique(summary_stats$scenario[i])),]$parental_contribution),
  R2<-lm(summary_stats[which(summary_stats$scenario==unique(summary_stats$scenario[i])),]$q.R2~summary_stats[which(summary_stats$scenario==unique(summary_stats$scenario[i])),]$parental_contribution)
  
}