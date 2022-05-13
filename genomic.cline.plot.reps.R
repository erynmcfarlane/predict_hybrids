###this is the first plot from composite.clines implemented in Introgress

genomic.cline.plot<-function(cline.data){
  plot(0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), cex.axis=1.5)
    for(i in 1:length(cline.data[[i]])){
      n.loci <- dim(cline.data[[i]]$Loci.data)[1]
      hi <- cline.data[[i]]$hybrid.index
      col<-c("red", "black", "black") ###hard coded for now, but could change depending on how many snps
      for (j in 1:n.loci) {
        AA.line <- cline.data[[i]]$Fitted.AA[j, ]
        line.matrix <- rbind(hi, AA.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 1, col=col[j])
        
        Aa.line <- cline.data[[i]]$Fitted.Aa[j, ]
        line.matrix <- rbind(hi, Aa.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 1, col=col[j])
  }
}
}
