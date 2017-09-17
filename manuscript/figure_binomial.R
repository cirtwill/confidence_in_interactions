library(binom)
n <- c(1:30)
CI <- binom.confint(0, n, method = "exact")
dev.new(width = 10, height = 8)
plot(n,CI[,6], xlab = "Number of samples", ylab = "Upper 95% confidence interval", type = "l",cex.axis = 1.5,cex.lab = 1.5,lwd = 2,ylim = c(0,1))
abline(h = 0.1, lty=3)


