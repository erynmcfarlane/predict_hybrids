#### Script for Figure 1 ####
source('Figure_prep.R')

triangle_cols <- met.brewer("OKeeffe1", 20)

quartz(height=8, width=12)
par(mar=c(5,5,0,0), oma=c(0,0,4,4), mfrow=c(4,6))

for (i in 1:24)
{
  if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme6_gen", anova_out_gen10[i,4], ".csv") }
  else							{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme6_gen", anova_out_gen10[i,4], ".csv") }
  dat_in <- read.csv(handle, header=FALSE)
  plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.75, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1, xaxt="n")
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col="light gray")
  axis(1, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=1.5)
  for (j in 1:20)
  {
    rep_sub <- subset(dat_in, dat_in[,1]==j)
    q_dens <- density(rep_sub[,5])
    points(q_dens$x, q_dens$y, type="l", lwd=1.5, col= triangle_cols[j])
  }
  if (i < 7)
  {
    mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
    mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
  }
  if (i %% 6 == 0)
  {
    if (anova_out_gen10[i,1]=="dmi") { mtext("DMI", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="dmi_e") { mtext("DMI + env", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="path") { mtext("Path", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="path_e") { mtext("Path + env", side=4, line=1.25, cex=1.5) }
  }
}
