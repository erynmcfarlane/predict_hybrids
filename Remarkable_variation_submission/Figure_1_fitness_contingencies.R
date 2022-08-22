#### Script for Figure 1 ####

## NOTE: dmis found at locations 1.4 and 1.6 (columns 12 & 14)


library(MetBrewer)
contin_cols <- met.brewer("Derain", 151)

contin <- read.csv("dmi_m0.01_c0.9_deme6_gen10.csv", header=F)
dim(contin)
contin[1:10,1:10]


## for 1.4 vs 1.6, columns should be set to 12 and 14
## for 2.4 vs 2.6, columns should be set to 63 and 65
## for 3.4 vs 3.6, columns should be set to 114 and 116

contin_mat <- matrix(0, 20, 9)
for (i in 1:20)
{
  rep_sub <- subset(contin, contin[,1]==i)
  for (j in 1:dim(rep_sub)[1])
  {
    if (rep_sub[j,12]==0)
    {
      if		(rep_sub[j,14]==0) { contin_mat[i,1] <- contin_mat[i,1] + 1 }
      else if	(rep_sub[j,14]==1) { contin_mat[i,2] <- contin_mat[i,2] + 1 }
      else if	(rep_sub[j,14]==2) { contin_mat[i,3] <- contin_mat[i,3] + 1 }
    }
    else if (rep_sub[j,12]==1)
    {
      if		(rep_sub[j,14]==0) { contin_mat[i,4] <- contin_mat[i,4] + 1 }
      else if	(rep_sub[j,14]==1) { contin_mat[i,5] <- contin_mat[i,5] + 1 }
      else if	(rep_sub[j,14]==2) { contin_mat[i,6] <- contin_mat[i,6] + 1 }
    }
    if (rep_sub[j,12]==2)
    {
      if		(rep_sub[j,14]==0) { contin_mat[i,7] <- contin_mat[i,7] + 1 }
      else if	(rep_sub[j,14]==1) { contin_mat[i,8] <- contin_mat[i,8] + 1 }
      else if	(rep_sub[j,14]==2) { contin_mat[i,9] <- contin_mat[i,9] + 1 }
    }
  }
}



## with fitness matrices

quartz(height=6, width=6)
layout(matrix(c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,
                1,1,1,1,1,1,1,2,2,2,2,2,2,2,
                1,1,1,1,1,1,1,2,2,2,2,2,2,2,
                3,3,4,4,5,5,6,6,7,7,8,8,9,9,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,
                17,17,18,18,19,19,20,20,21,21,22,22,23,23), 6, 14, byrow=TRUE))
par(mar=c(3,3,2,2))

## DMI fitness matrix
plot(0, type="n", axes=FALSE, xlab="", ylab="", main="DMI", xlim=c(0,3), ylim=c(0,3), xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5)
axis(2, at=c(2.5, 1.5, 0.5), labels=c("aa", "aA", "AA"), cex.axis=1.5, lwd=0, line=-1, las=1)
axis(1, at=c(2.5, 1.5, 0.5), labels=c("bb", "bB", "BB"), cex.axis=1.5, lwd=0, line=-1)
rect(0, 0, 1, 1, lwd=2)
rect(1, 0, 2, 1, lwd=2)
rect(2, 0, 3, 1, lwd=2)
rect(0, 1, 1, 2, lwd=2)
rect(1, 1, 2, 2, lwd=2)
rect(2, 1, 3, 2, lwd=2)
rect(0, 2, 1, 3, lwd=2)
rect(1, 2, 2, 3, lwd=2)
rect(2, 2, 3, 3, lwd=2)
text(0.5, 2.5, labels="1", cex=1.5)
text(0.5, 1.5, labels="1 - c", cex=1.5)
text(0.5, 0.5, labels="1 - c", cex=1.5)
text(1.5, 2.5, labels="1", cex=1.5)
text(1.5, 1.5, labels="1 - c", cex=1.5)
text(1.5, 0.5, labels="1 - c", cex=1.5)
text(2.5, 2.5, labels="1", cex=1.5)
text(2.5, 1.5, labels="1", cex=1.5)
text(2.5, 0.5, labels="1", cex=1.5)
mtext("A", adj=-0.15, cex=1.5)

## Pathways fitness matrix
plot(0, type="n", axes=FALSE, xlab="", ylab="", main="Pathway", xlim=c(0,3), ylim=c(0,3), xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5)
axis(2, at=c(2.5, 1.5, 0.5), labels=c(expression("A"[1]~"A"[1]), expression("A"[1]~"A"[2]), expression("A"[2]~"A"[2])), cex.axis=1.5, lwd=0, line=-1, las=1)
axis(1, at=c(2.5, 1.5, 0.5), labels=c(expression("B"[2]~"B"[2]), expression("B"[1]~"B"[2]), expression("B"[1]~"B"[1])), cex.axis=1.5, lwd=0, line=-1)
rect(0, 0, 1, 1, lwd=2)
rect(1, 0, 2, 1, lwd=2)
rect(2, 0, 3, 1, lwd=2)
rect(0, 1, 1, 2, lwd=2)
rect(1, 1, 2, 2, lwd=2)
rect(2, 1, 3, 2, lwd=2)
rect(0, 2, 1, 3, lwd=2)
rect(1, 2, 2, 3, lwd=2)
rect(2, 2, 3, 3, lwd=2)
text(0.5, 2.5, labels="1", cex=1.5)
text(0.5, 1.5, labels="1 - c/2", cex=1.5)
text(0.5, 0.5, labels="1 - c", cex=1.5)
text(1.5, 2.5, labels="1 - c/2", cex=1.5)
text(1.5, 1.5, labels="1", cex=1.5)
text(1.5, 0.5, labels="1 - c/2", cex=1.5)
text(2.5, 2.5, labels="1 - c", cex=1.5)
text(2.5, 1.5, labels="1 - c/2", cex=1.5)
text(2.5, 0.5, labels="1", cex=1.5)
mtext("B", adj=-0.15, cex=1.5)

## contingency
par(mar=c(1, 0.5, 1, 0.5))
for (i in 1:20)
{
  plot(0, type="n", xlim=c(0,3), ylim=c(0,3), xaxt="n", yaxt="n", xlab="", ylab="", main=paste0("rep ", i), axes=FALSE)
  #axis(2, at=c(0.5, 1.5, 2.5), labels=c(2, 1, 0), las=1, cex.axis=2)
  #axis(3, at=c(0.5, 1.5, 2.5), labels=c(0, 1, 2), cex.axis=2)
  rect(0, 2, 1, 3, col=contin_cols[(contin_mat[i,1] + 1)], lwd=2)
  rect(1, 2, 2, 3, col=contin_cols[(contin_mat[i,2] + 1)], lwd=2)
  rect(2, 2, 3, 3, col=contin_cols[(contin_mat[i,3] + 1)], lwd=2)
  rect(0, 1, 1, 2, col=contin_cols[(contin_mat[i,4] + 1)], lwd=2)
  rect(1, 1, 2, 2, col=contin_cols[(contin_mat[i,5] + 1)], lwd=2)
  rect(2, 1, 3, 2, col=contin_cols[(contin_mat[i,6] + 1)], lwd=2)
  rect(0, 0, 1, 1, col=contin_cols[(contin_mat[i,7] + 1)], lwd=2)
  rect(1, 0, 2, 1, col=contin_cols[(contin_mat[i,8] + 1)], lwd=2)
  rect(2, 0, 3, 1, col=contin_cols[(contin_mat[i,9] + 1)], lwd=2)
  if (i==1) { mtext("C", adj=0, cex=1.5, line=1) }
}

## legend
par(mar=c(1, 1, 1, 1))
legend_image <- as.raster(matrix(rev(contin_cols), ncol=1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="", main="# individuals")
rasterImage(legend_image, 1/6, 0.05, 1/2, 0.95)
rect(1/6, 0.05, 1/2, 0.95, lwd=1.25)
text(0.75, 0.1, labels="0", cex=1.25)
text(0.75, 0.5, labels="75", cex=1.25)
text(0.75, 0.9, labels="150", cex=1.25)


