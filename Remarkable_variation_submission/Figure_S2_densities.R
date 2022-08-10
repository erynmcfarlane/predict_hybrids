#### Script for Figure S2 ###
source('Figure_prep.R')
source('hexagon.R')

quartz(height=8, width=20)
layout(matrix(c(25,25,25,25,1,2,3,4,5,6,
                25,25,25,25,7,8,9,10,11,12,
                25,25,25,25,13,14,15,16,17,18,
                25,25,25,25,19,20,21,22,23,24), 4, 10, byrow=TRUE))
par(mar=c(5,5,0,0), oma=c(0,0,4,4))

for (i in 1:24)
{
  if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme6_gen", anova_out_gen10[i,4], ".csv") }
  else										{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme6_gen", anova_out_gen10[i,4], ".csv") }
  dat_in <- read.csv(handle, header=FALSE)
  plot(0, type="n", xlab="q", ylab="Density", cex.lab=1.5, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1)
  if		(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.2)		{ points(0.9, 10.8, pch=21, bg=anova_cols[i], cex=3) }
  else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.9)		{ points(0.9, 10.8, pch=22, bg=anova_cols[i], cex=3) }
  #else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]=="neutral")	{ points(0.9, 10.8, pch=13, col=anova_cols[i], cex=3) }
  else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0)	{ hexagon(x_center=0.9, y_center=10.8, scale=0.075, plot_xmin=0, plot_xmax=1, plot_ymin=0, plot_ymax=12, color=anova_cols[i]) }
  else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.2)			{ points(0.9, 10.8, pch=24, bg=anova_cols[i], cex=3) }
  else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.9)			{ points(0.9, 10.8, pch=25, bg=anova_cols[i], cex=3) }
  else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0)	{ points(0.9, 10.8, pch=23, bg=anova_cols[i], cex=3) }
  for (j in 1:20)
  {
    rep_sub <- subset(dat_in, dat_in[,1]==j)
    q_dens <- density(rep_sub[,5])
    points(q_dens$x, q_dens$y, type="l", lwd=0.5, col=anova_cols[i])
    polygon(q_dens, col=adjustcolor(anova_cols[i], alpha.f=0.2), border=FALSE)
  }
  if (i < 7)
  {
    mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
    mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
  }
  if (i %% 6 == 0)
  {
    if (anova_out_gen10[i,1]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
    else if (anova_out_gen10[i,1]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
  }
}
plot(anova_out_gen10[,5], anova_out_gen10[,6], type="n", xlab="q anova log10 F", ylab="Q anova log10 F", cex.lab=2.5, cex.axis=1.5, las=1, xlim=range(as.numeric(anova_out_gen10[,5:6])), ylim=range(as.numeric(anova_out_gen10[,5:6])))
abline(0, 1, lty=2, lwd=4)
for (i in 1:24)
{
  if		(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.2)		{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=21, bg=anova_cols[i], cex=5) }
  else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0.9)		{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=22, bg=anova_cols[i], cex=4) }
  else if	(anova_out_gen10[i,2]==0.01 & anova_out_gen10[i,3]==0)	{ hexagon(as.numeric(anova_out_gen10[i,5]), as.numeric(anova_out_gen10[i,6]), scale=0.02, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color=anova_cols[i]) }
  else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.2)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=24, bg=anova_cols[i], cex=4) }
  else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0.9)			{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=25, bg=anova_cols[i], cex=4) }
  else if	(anova_out_gen10[i,2]==0.2 & anova_out_gen10[i,3]==0)	{ points(anova_out_gen10[i,5], anova_out_gen10[i,6], pch=23, bg=anova_cols[i], cex=4) }
}
legend("topleft", legend=c("m = 0.01; c = 0", "m = 0.01; c = 0.2", "m = 0.01; c = 0.9", "m = 0.2; c = 0", "m = 0.2; c = 0.2", "m = 0.2; c = 0.9"), pch=c(20,21,22,23,24,25), pt.cex=3, pt.bg="black", cex=2)
hexagon(-0.178, 3.06, scale=0.015, plot_xmin=min(as.numeric(anova_out_gen10[,5:6])), plot_xmax=max(as.numeric(anova_out_gen10[,5:6])), plot_ymin=min(as.numeric(anova_out_gen10[,5:6])), plot_ymax=max(as.numeric(anova_out_gen10[,5:6])), color="black")




