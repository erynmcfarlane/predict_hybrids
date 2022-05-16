###this is the first plot from composite.clines implemented in Introgress

genomic.cline.plot<-function(cline.data){
    for(k in 1:length(cline.data)){
      n.loci <- dim(cline.data[[k]]$Loci.data)[1]
      hi <- cline.data[[k]]$hybrid.index
      col<-c("red", "black", "black") ###hard coded for now, but could change depending on how many snps
      for (j in 1:n.loci) {
        AA.line <- cline.data[[k]]$Fitted.AA[j, ]
        line.matrix <- rbind(hi, AA.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 1, col=col[j])
        
        Aa.line <- cline.data[[k]]$Fitted.Aa[j, ]
        line.matrix <- rbind(hi, Aa.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 1, col=col[j])
  }
}
}
