###this is the first plot from composite.clines implemented in Introgress

genomic.cline.plot<-function(cline.data, snp.no){ ### I've now put this in to tell it which SNP i want it to plot so that I don't have to have it hard coded in.
    for(k in 1:length(cline.data)){
      n.loci <- dim(cline.data[[k]]$Loci.data)[1]
      hi <- cline.data[[k]]$hybrid.index
      col<-colours 

      for (l in snp.no:snp.no) {
         AA.line <- cline.data[[k]]$Fitted.AA[l, ]
        line.matrix <- rbind(hi, AA.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 1,lwd=4, col=col[k])
        
        Aa.line <- cline.data[[k]]$Fitted.Aa[l, ]
        line.matrix <- rbind(hi, Aa.line)[, order(hi)]
        lines(line.matrix[1, ], line.matrix[2, ], lty = 2, lwd=4,col=col[k])
  }
}
}
