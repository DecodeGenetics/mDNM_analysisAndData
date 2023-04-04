mutRateData <- read.table('./mutRateDataAll', header = TRUE)
glm_poiss_offs <- glm(DNMs ~ mutRateData$motifLength + refLen + purity + GCcontentInMotif + offset(log(Correct)), family=poisson(link='log'), data=mutRateData)

motifLength <- c(1,2,3,4,5,6)
pvals <- matrix(nrow = 6, ncol = 4)
coefs <- matrix(nrow = 6, ncol = 4)
for (i in motifLength) {
  glm_poiss_offs <- glm(DNMs ~ refLen + purity + GCcontentInMotif + offset(log(Correct)), family=poisson(link='log'), data=mutRateData[mutRateData$motifLength==i,])
  pvals[i,] <- coef(summary(glm_poiss_offs))[,4]
  coefs[i,] <- glm_poiss_offs$coefficients
}

refLens <- c(10,20,30,40,50,60,70,80,90)
pvals_2 <- matrix(nrow = 9, ncol = 6)
coefs_2 <- matrix(nrow = 9, ncol = 6)
for (i in motifLength) {
  for (j in refLens)
  {
    if (dim(mutRateData[mutRateData$motifLength==i & mutRateData$refLen>j & mutRateData$refLen<(j+10),])>1)
    {
      glm_poiss_offs <- glm(DNMs ~ purity + offset(log(Correct)), family=poisson(link='log'), data=mutRateData[mutRateData$motifLength==i & mutRateData$refLen>j & mutRateData$refLen<(j+10),])
      pvals_2[j/10,i] <- coef(summary(glm_poiss_offs))[,4][2]
      coefs_2[j/10,i] <- glm_poiss_offs$coefficients[2]
    }
    else
    {
      pvals_2[j/10,i] <- NA
      coefs_2[j/10,i] <- NA
    }
  }
}