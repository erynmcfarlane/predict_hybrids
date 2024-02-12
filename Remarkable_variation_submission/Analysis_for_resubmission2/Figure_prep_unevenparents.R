###Prep script for most figures ###
library(MetBrewer)
anova_cols <- met.brewer("Hokusai1", 12)

uniq_runs <- c("dmi_m", "path_m")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)
uniq_gen <- c(10, 100)

anova_out <- matrix(0, 24, 13)

ctr <- 1
for (i in 1:length(uniq_runs))
{
  for (j in 1:length(uniq_m))
  {
    for (k in 1:length(uniq_c))
    {
      for (l in 1:length(uniq_gen))
      {
        print(ctr)
        dat_in <- data_deme_2[which(data_deme_2$mech==uniq_runs[i] & data_deme_2$m==uniq_m[j] & data_deme_2$c==uniq_c[k] & data_deme_2$gen==uniq_gen[l]),]
        q_anova <- aov(dat_in$q ~ as.factor(dat_in$rep))
        het_anova <- aov(dat_in$het ~ as.factor(dat_in$rep))
        junc_anova <- aov(dat_in$njunct ~ as.factor(dat_in$rep))
        
        anova_out[ctr,1] <- uniq_runs[i]
        anova_out[ctr,2] <- uniq_m[j]
        anova_out[ctr,3] <- uniq_c[k]
        anova_out[ctr,4] <- uniq_gen[l]
        anova_out[ctr,5] <- log10(summary(q_anova)[[1]][1,4])
        anova_out[ctr,6] <- log10(summary(het_anova)[[1]][1,4])
        anova_out[ctr,7] <- log10(summary(junc_anova)[[1]][1,4])
        anova_out[ctr,8] <- summary(q_anova)[[1]][1,4]
        anova_out[ctr,9] <- summary(q_anova)[[1]][1,5]
        anova_out[ctr,10] <- summary(het_anova)[[1]][1,4]
        anova_out[ctr,11] <- summary(het_anova)[[1]][1,5]
        anova_out[ctr,12] <- summary(junc_anova)[[1]][1,4]
        anova_out[ctr,13] <- summary(junc_anova)[[1]][1,5]
        ctr <- ctr+1
      }
    }
  }
}

## export for table
table_anova_out <- anova_out[,c(1:4,8:13)]
colnames(table_anova_out) <- c("run", "m", "c", "gen", "q_f", "q_p", "het_f", "het_p", "junc_f", "junc_p")
#write.table(table_anova_out, file="anova_out_table.txt", row.names=F, quote=F, sep="\t") ### what is the table that is actually exported called? There should be 12 tables?

anova_out_gen10 <- subset(anova_out, anova_out[,4]=="10")
anova_out_gen100 <- subset(anova_out, anova_out[,4]=="100")