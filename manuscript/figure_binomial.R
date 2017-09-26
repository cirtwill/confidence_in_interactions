library(binom)
n <- c(1:30)
CI <- binom.confint(0, n, method = "exact")
dev.new(width = 10, height = 8)
plot(n,CI[,6], xlab = "Number of samples", ylab = "Upper 95% confidence interval", type = "l",cex.axis = 1.5,cex.lab = 1.5,lwd = 2,ylim = c(0,1))
abline(h = 0.1, lty=3)
write.table(CI,file='../data/R/CI_and_samples.tsv',sep='\t')
# CAPTION
# This figure illustrates the upper confidence interval for the interaction probability $\lamba$ when there is no observation of an interaction ($k = 0$) and $n$ observations of the two species co-occurring at a given location. The confidence interval is based on the Clopper-Pearson derivation from the beta distribution. 


