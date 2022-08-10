library("vegan")
uniq_runs <- c("dmi", "dmi_e", "path", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)

gen <- 100
quartz(height=12, width=18)
par(mar=c(5,5,0,0), mfrow=c(4,6), oma=c(0,0,4,4))
ctr <- 0
for (i in 1:length(uniq_runs))
{
  for (j in 1:length(uniq_m))
  {
    for (k in 1:length(uniq_c))
    {
      ctr <- ctr + 1
      if (uniq_c[k]==0)	{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen", gen, ".csv") }
      else					{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen", gen, ".csv") }
      dat_in <- read.csv(handle, header=FALSE)
      pca_out <- prcomp(dat_in[,9:518], center=TRUE, scale=FALSE)
      pc1_exp <- round(summary(pca_out)$importance[2,1]*100, 2)
      pc2_exp <- round(summary(pca_out)$importance[2,2]*100, 2)
      plot(pca_out$x[,1], pca_out$x[,2], type="n", xlab=paste0("PC1 (", pc1_exp, "%)"), ylab=paste0("PC2 (", pc2_exp, "%)"), cex.lab=1.75, cex.axis=1.5, las=1, xlim=c(min(pca_out$x[,1]*2), max(pca_out$x[,1]*2)), ylim=c(min(pca_out$x[,2]*1.1), max(pca_out$x[,2]*1.1))); box(lwd=2)
      pca_plot_data <- cbind(dat_in[,1:4], pca_out$x[,1:2])
      points(pca_out$x[,1], pca_out$x[,2], col="gray", cex=0.5)
      ordiellipse(pca_out, pca_plot_data[,1], conf=0.9, col="red")
      if (ctr < 7)
      {
        mtext(paste0("m = ", uniq_m[j]), line=2, cex=1.5)
        mtext(paste0("c = ", uniq_c[k]), line=0.25, cex=1.5)
      }
      if (ctr %% 6 == 0)
      {
        if 		(uniq_runs[i]=="dmi") { mtext("dmi", side=4, line=1.25, cex=1.5) }
        else if (uniq_runs[i]=="dmi_e") { mtext("dmi + env", side=4, line=1.25, cex=1.5) }
        else if (uniq_runs[i]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
        else if (uniq_runs[i]=="path_e") { mtext("path + env", side=4, line=1.25, cex=1.5) }
      }
    }
  }
}