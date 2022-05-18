###this is the first plot from composite.clines implemented in Introgress

genomic.cline.plot<-function(cline.data){
    for(k in 1:length(cline.data)){
      n.loci <- dim(cline.data[[k]]$Loci.data)[1]
      hi <- cline.data[[k]]$hybrid.index
      col<-colours 
      ###hard coded for now, but could change depending on how many snps
      for (l in 1:length(n.loci)) {
         AA.line <- cline.data[[k]]$Fitted.AA[l, ]
        line.matrix <- rbind(hi, AA.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 1, col=col[k])
        
        Aa.line <- cline.data[[k]]$Fitted.Aa[l, ]
        line.matrix <- rbind(hi, Aa.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 2, col=col[k])
  }
}
}
