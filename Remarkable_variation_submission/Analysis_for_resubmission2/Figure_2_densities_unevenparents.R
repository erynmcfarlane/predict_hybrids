#### Script for Figure 2 ####
source('SNP_inputs_unevenparents.R')
source('Figure_prep_unevenparents.R')
### instead of anova_out_gen10, this is alldata_df

triangle_cols <- met.brewer("OKeeffe1", 20)

#quartz(height=4, width=12)
pdf(height=4, width=12, file="anova_gen10_rep_lines_twoone.pdf")
par(mar=c(6,6,0,0), oma=c(0,0,4,4), mfrow=c(2,6))

###anova_out_gen10 has only 12 rows in it. Need to read in the actualy data 1:12 based on this. 
### so need to paste together the categories from anova_out_gen10 to pull from data_deme_2
##data_deme_2[which(data_deme_2$mech==anova_out_gen10$[i,1] & data_deme_2$m==anova_out_gen10$[i,2] & data_deme_2$c==anova_out_gen10$[i,3] & data_deme_2$gen==anova_out_gen10$[i,4],]


for (i in 1:12) ### clearly, this is still not working.
{
  dat_in <- data_deme_2[which(data_deme_2$mech==anova_out_gen10[i,1] & data_deme_2$m==anova_out_gen10[i,2] & data_deme_2$c==anova_out_gen10[i,3] & data_deme_2$gen==anova_out_gen10[i,4]), ]
  plot(0, type="n", xlab="", ylab="", cex.lab=1.75, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1, xaxt="n")
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col="light gray")
  axis(1, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=1.5)
  for (j in 1:20)
  {
    rep_sub <- dat_in[which(dat_in$rep==j),]
    q_dens <- density(rep_sub$q)
    points(q_dens$x, q_dens$y, type="l", lwd=1.5, col= triangle_cols[j])
  }
  if (i < 7)
  {
    mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
    mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
  }
  if (i %% 6 == 0)
  {
    if (anova_out_gen10[i,1]=="dmi") { mtext("BDMI", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
  }
  if (i==10) { mtext("Admixture proportion (q)", side=1, line=4, adj=0.95, cex=2) }
  if (i==7) { mtext("Density", side=2, line=3.5, adj=-5, cex=2) }
}

dev.off()




