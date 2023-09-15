########################################
## 1) load data and packages
########################################

library(MetBrewer)




########################################
## 11) anova plots
########################################


uniq_runs <- c("dmi", "path")
uniq_m <- c(0.01, 0.2)
uniq_c <- c(0, 0.2, 0.9)
uniq_gen <- c(10)

anova_out <- matrix(0, 12, 13)

ctr <- 1
for (i in 1:length(uniq_runs))
	{
	for (j in 1:length(uniq_m))
		{
		for (k in 1:length(uniq_c))
			{
			for (l in 1:length(uniq_gen))
				{
				print(ctr)
				if (uniq_c[k]==0)	{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_neutral", "_deme2_gen", uniq_gen[l], ".csv") }
				else							{ handle <- paste0(uniq_runs[i], "_m", uniq_m[j], "_c", uniq_c[k], "_deme2_gen", uniq_gen[l], ".csv") }
				dat_in <- read.csv(handle, header=FALSE)
				q_anova <- aov(dat_in[,5] ~ as.factor(dat_in[,1]))
				het_anova <- aov(dat_in[,6] ~ as.factor(dat_in[,1]))
				junc_anova <- aov(dat_in[,8] ~ as.factor(dat_in[,1]))
				
				anova_out[ctr,1] <- uniq_runs[i]
				anova_out[ctr,2] <- uniq_m[j]
				anova_out[ctr,3] <- uniq_c[k]
				anova_out[ctr,4] <- uniq_gen[l]
				anova_out[ctr,5] <- log10(summary(q_anova)[[1]][1,4])
				anova_out[ctr,6] <- log10(summary(het_anova)[[1]][1,4])
				anova_out[ctr,7] <- log10(summary(junc_anova)[[1]][1,4])
				anova_out[ctr,8] <- summary(q_anova)[[1]][1,4]
				anova_out[ctr,9] <- summary(q_anova)[[1]][1,5]
				anova_out[ctr,10] <- summary(het_anova)[[1]][1,4]
				anova_out[ctr,11] <- summary(het_anova)[[1]][1,5]
				anova_out[ctr,12] <- summary(junc_anova)[[1]][1,4]
				anova_out[ctr,13] <- summary(junc_anova)[[1]][1,5]
				ctr <- ctr+1
				}
			}
		}
	}
	
## export for table
table_anova_out <- anova_out[,c(1:4,8:13)]
colnames(table_anova_out) <- c("run", "m", "c", "gen", "q_f", "q_p", "het_f", "het_p", "junc_f", "junc_p")
#write.table(table_anova_out, file="anova_out_table.txt", row.names=F, quote=F, sep="\t")

anova_out_gen10 <- subset(anova_out, anova_out[,4]=="10")



## test different colors for different replicates

triangle_cols <- met.brewer("OKeeffe1", 20)

quartz(height=4, width=12)
par(mar=c(3,3,0,0), oma=c(4,4,4,4), mfrow=c(2,6))

for (i in 1:12)
	{
	if (anova_out_gen10[i,3]==0)	{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_neutral", "_deme2_gen", anova_out_gen10[i,4], ".csv") }
	else							{ handle <- paste0(anova_out_gen10[i,1], "_m", anova_out_gen10[i,2], "_c", anova_out_gen10[i,3], "_deme2_gen", anova_out_gen10[i,4], ".csv") }
	dat_in <- read.csv(handle, header=FALSE)
	plot(0, type="n", xlab="", ylab="", cex.lab=1.75, cex.axis=1.5, ylim=c(0,12), xlim=c(0,1), las=1, xaxt="n")
	rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col="light gray")
	axis(1, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=1.5)
	for (j in 1:20)
		{
		rep_sub <- subset(dat_in, dat_in[,1]==j)
		q_dens <- density(rep_sub[,5])
		points(q_dens$x, q_dens$y, type="l", lwd=1.5, col= triangle_cols[j])
		}
	if (i < 7)
		{
		mtext(paste0("m = ", anova_out_gen10[i,2]), line=2, cex=1.5)
		mtext(paste0("c = ", anova_out_gen10[i,3]), line=0.25, cex=1.5)
		}
	if (i %% 6 == 0)
		{
		if (anova_out_gen10[i,1]=="dmi") { mtext("BDMI", side=4, line=1.25, cex=1.5) }
		else if (anova_out_gen10[i,1]=="path") { mtext("path", side=4, line=1.25, cex=1.5) }
		}
	if (i==10) { mtext("Admixture proportion (q)", side=1, line=4, adj=0.95, cex=2) }
	if (i==1) { mtext("Density", side=2, line=4, adj=22, cex=2) }
	}
