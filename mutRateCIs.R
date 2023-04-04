library(boot)

motifLength <- c(1,2,3,4,5,6)

#function to compute mutation rate
function_1 <- function(data, i){
  d2 <- data[i,]
  return(sum(d2[,6])/(sum(d2[,6])+sum(d2[,7])))
}

mutRateData <- read.table('./mutRateDataAll', header = TRUE)
mutRateData$motif <- as.character(mutRateData$motif)
bootstrap_mutRate <- boot(mutRateData, function_1, R=100)
#vectors to store mutation rate and confidence intervals on motif length
mutRate_motifLength <- rep(NA, 6)
mutRate_motifLength_LL <- rep(NA, 6)
mutRate_motifLength_UL <- rep(NA, 6) 
#mut rate per motif length with CI
for (i in motifLength) {
  bootstrap_motifLength <- boot(mutRateData[nchar(mutRateData$motif) == i, ], function_1, R=100)
  mutRate_motifLength[i] <- bootstrap_motifLength$t0
  mutRate_motifLength_LL[i] <- quantile(bootstrap_motifLength$t, probs=seq(0,1,0.025))[[2]]
  mutRate_motifLength_UL[i] <- quantile(bootstrap_motifLength$t, probs=seq(0,1,0.025))[[40]]
}