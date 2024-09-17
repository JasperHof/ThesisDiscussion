
library(ggplot2)

### create sets p1, p2, p3 (calibrated, inflated (relatedness), inflated (spa))

n_snps = 5e5

### 1

p1 = runif(n_snps)

### 2

N = 750
G = matrix(rbinom(N * n_snps, 2, 0.05), nrow = N)
G = rbind(G, G[1:250,])

preds = sample(1:n_snps, 1000)
Y = G[,preds] %*% rnorm(1000); Y = as.numeric(scale(Y))

lower = apply(G, 2, function(x) sum(x^2))
upper = apply(G, 2, function(x) sum(x*Y))

beta = upper/lower
beta_var = (sum(Y^2)/1000)/(lower)
stats = beta/sqrt(beta_var)

p2 = 2*(1 - pnorm(abs(stats)))

### 3

N = 1000
G = matrix(rbinom(N * n_snps, 2, 0.05), nrow = N)

preds = sample(1:n_snps, 1000)
liab = G[,preds] %*% rnorm(1000); liab = as.numeric(scale(liab))
Y = rbinom(N, 1, exp(-liab)/(1+exp(3-liab))); Y = scale(Y)

lower = apply(G, 2, function(x) sum(x^2))
upper = apply(G, 2, function(x) sum(x*Y))

beta = upper/lower
beta_var = (sum(Y^2)/1000)/(lower)
stats = beta/sqrt(beta_var)

p3 = 2*(1 - pnorm(abs(stats)))

### plot

all_p = cbind(p1[order(p1)], p2[order(p2)], p3[order(p3)], ppoints(length(p1)))
all_p = as.data.frame(all_p)
colnames(all_p) = c('p1','p2','p3','qq')

# reduce this a bit
all_p = all_p[-sample(which(all_p$qq > 0.05), floor(sum(all_p$qq > 0.05)*9.5/10)),]
all_p = all_p[-sample(which(all_p$qq > 0.001), floor(sum(all_p$qq > 0.001)*8/10)),]
all_p = apply(all_p,2, function(x) -log10(x))
all_p = as.data.frame(all_p)
colnames(all_p) = c('p1','p2','p3','qq')

par(mfrow = c(1,3))

plot(all_p$qq, all_p$p1, ylim=c(0,8), main = 'No inflation', xlab = 'Expected quantiles', ylab = 'Observed quantiles')
abline(0,1, col = 'red')
text(0.5,7,expression(bold('A')), cex  = 1.5)
plot(all_p$qq, all_p$p3, ylim=c(0,8), main = 'Inflation due to small sample size\n or case-control imbalance', xlab = 'Expected quantiles', ylab = 'Observed quantiles')
abline(0,1, col = 'red')
text(0.5,7,expression(bold('B')), cex  = 1.5)
plot(all_p$qq, all_p$p2, ylim=c(0,8), main = 'Inflation due to relatedness', xlab = 'Expected quantiles', ylab = 'Observed quantiles')
abline(0,1, col = 'red')
text(0.5,7,expression(bold('C')), cex  = 1.5)

save.image('inflation_p.R')

################### CONTINUE FROM HERE #######################

load('inflation_p.R')

library(ggplot2)

# determine confidence intervals for all curves
head(all_p)
#all_p$lower = 

p1 = ggplot(all_p, aes(qq, p1)) + ylim(0,8) + geom_abline(intercept = 0, slope = 1, col = 'red') + geom_point() +
   xlab('Expected quantiles') + ylab('Observed quantiles') + ggtitle('No inflation') + theme_light() + 
   theme(
     plot.title = element_text(size  = 20,hjust = 0.5),
     axis.title = element_text(size = 18),
     axis.text = element_text(size = 16)
   ) 
p2 = ggplot(all_p, aes(qq, p3)) + ylim(0,8) + geom_abline(intercept = 0, slope = 1, col = 'red') + geom_point() +
  xlab('Expected quantiles') + ylab('') + ggtitle('Small sample size/case-control imbalance') + theme_light() + 
  theme(
    plot.title = element_text(size  = 20,hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) 
p3 = ggplot(all_p, aes(qq, p2)) + ylim(0,8) + geom_abline(intercept = 0, slope = 1, col = 'red') + geom_point() +
  xlab('Expected quantiles') + ylab('') + ggtitle('Relatedness') + theme_light() + 
  theme(
    plot.title = element_text(size  = 20,hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) 
  
library(gridExtra)
  
png("//umcsanfsclp01/hev_genepi/genepi/Jasper/Papers0_Thesis/Figures/inflation.png",
    width = 15, height = 5, units = 'in', res = 1200)

grid.arrange(p1,p2,p3, nrow=1)

dev.off()

### add confidence intervals

