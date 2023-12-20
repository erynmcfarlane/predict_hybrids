### Script for Figure 4 ###

library(MetBrewer)
triangle_cols <- met.brewer("OKeeffe1", 20)

triangle <- read.csv("dmi_m0.01_c0.9.main_first8cols.txt", header=TRUE)
dim(triangle)
head(triangle)

gen10 <- subset(triangle, triangle[,2]==10)
dim(gen10)
head(gen10)

gen10_deme6 <- subset(gen10, gen10[,3]==6)
dim(gen10_deme6)
head(gen10_deme6)

quartz(height=12, width=15)
par(mar=c(5,5,1,1), mfrow=c(4,5))
for (i in 1:20)
{
  rep_sub <- subset(gen10_deme6, gen10_deme6[,1]==i)
  plot(0, type="n", xlab="q", ylab="Q", xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
  segments(0, 0, 0.5, 1, lwd=3)
  segments(1, 0, 0.5, 1, lwd=3)
  for (j in 1:150)
  {
    points(rep_sub[j,5], rep_sub[j,6], pch=21, bg=adjustcolor(triangle_cols[i], alpha.f=0.75), cex=2)
  }
  box(lwd=2)
}