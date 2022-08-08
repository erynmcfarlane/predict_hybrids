###this is the first plot from composite.clines implemented in Introgress
### I either need to get Josh's help to make these
genomic.cline.plot<-function(cline.data){
  #plot(0:1, 0:1, type = "n", xlab = "Hybrid index", ylab = 'Probability of Ancestry',cex.lab=2, main=main_label )
 cline.data<-as.data.frame(cline.data)
  n.loci <- dim(cline.data$Loci.data)[1]
  hi <- cline.data$hybrid.index
  col<-c("red", "black", "black") ###hard coded for now, but could change depending on how many snps
  for (i in 1:n.loci) {
    AA.line <- cline.data$Fitted.AA[i, ]
    line.matrix <- rbind(hi, AA.line)[, order(hi)]
    ggplot(cline.data, aes(x=hybrid.index, y=Fitted.AA[i,]))+geom_line(line.matrix[1, ], line.matrix[2, ])
  }
}
cline.data<-as.data.frame(genomic.clines[[i]])


ahhhhhhhh
