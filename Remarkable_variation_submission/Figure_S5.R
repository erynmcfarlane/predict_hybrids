### Script for Figure S5 ###
library(MetBrewer)
rep_cols <- met.brewer("OKeeffe1", 20)


uniq_runs <- c("dmi", "dmi_e", "path", "path_e")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)

quartz(height=12, width=18)
par(mar=c(5,5,1,1), mfrow=c(4,6), oma=c(0,0,4,4))
ctr <- 0
for (i in 1:length(uniq_runs))
{
  for (j in 1:length(uniq_m))
  {
    for (k in 1:length(uniq_c))
    {
      ctr <- ctr + 1
      if (uniq_c[k]==0)
      {
        handle10 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen10.csv")
        handle100 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme6_gen100.csv")
      }
      else	
      {
        handle10 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen10.csv")
        handle100 <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme6_gen100.csv")
      }
      
      dat10 <- read.csv(handle10, header=FALSE)
      dat100 <- read.csv(handle100, header=FALSE)
      plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="Median gen 10 Q", ylab="Median gen 100 Q", cex.lab=1.75, cex.axis=1.5, las=1)
      box(lwd=2)
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
      cor_mat <- matrix(0, 20, 2)
      for (l in 1:20)
      {
        rep_sub10 <- subset(dat10, dat10[,1]==l)
        rep_sub100 <- subset(dat100, dat100[,1]==l)
        cor_mat[l,1] <- median(rep_sub10[,6])
        cor_mat[l,2] <- median(rep_sub100[,6])
        points(median(rep_sub10[,6]), median(rep_sub100[,6]), pch=21, bg=rep_cols[l], cex=2)
      }
      mtext(paste0("r = ", round(cor(cor_mat[,1], cor_mat[,2], method="pearson"), 3)), line=-2, adj=0.05, cex=1.25)
    }
  }
}