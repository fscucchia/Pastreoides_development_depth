
#Visualize the extracted information, plot DP distribution and define cut-off
#See original script here https://evodify.com/gatk-in-non-model-organism/

# define the file names list 
setwd("/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration")


nameList <- c()
for (i in 3:25) { # 25 - odd number for 12 samples 
  if (i %% 2 ==1) nameList <- append(nameList, paste(i, ".DP"))
}

qlist <- matrix(nrow = 12, ncol = 3) # define number of samples (10 samples here)
qlist <- data.frame(qlist, row.names=nameList)
colnames(qlist)<-c('5%', '10%', '99%')

jpeg("GVCFall.DP.jpeg", height=1600, width=1200)
par(mar=c(5, 3, 3, 2), cex=1.5, mfrow=c(10,6)) # define number of plots for your sample
for (i in 1:37) {
  DP <- read.table(nameList[i], header = T)
  qlist[i,] <- quantile(DP[,1], c(.05, .1, .99), na.rm=T)
  d <- density(DP[,1], from=0, to=100, bw=1, na.rm =T)
  plot(d, xlim = c(0,100), main=nameList[i], col="blue", xlab = dim(DP)[1], lwd=2)
  abline(v=qlist[i,c(1,3)], col='red', lwd=3)
}
dev.off()

write.table(qlist, "GVCFall.DP.percentiles.txt")