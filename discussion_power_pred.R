
library(readxl)
library(survminer)
library(survival)
library(simrec)
library(tidyr)

### Figures:

# distribution                         bin, tte, tr, quant
# power (mean chisq, effs 0 - 1)       bin, tte, tr, quant
# prs   

all = NULL

### introduction - sample sizes?

for(j in 1:100){

N = 1000
n_snps = 22

G = matrix(rbinom(N * n_snps, 2, 0.5), nrow = N)
eff = rep(seq(from = 0, to = 1, by = 0.1), each = 2); eff = c(eff, rep(0,n_snps - length(eff)))
g_liab = G %*% eff;

# introduce some frailty (for RE)
frac = 0.8
frail = rnorm(N)

liab = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab = scale(liab)    # This is already the quantitative phenotype
liab2 = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab2 = scale(liab2) #2nd recurrence
liab3 = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab3 = scale(liab3) #3rd recurrence
liab4 = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab4 = scale(liab4) #4th recurrence
liab5 = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab5 = scale(liab5) #5nd recurrence
liab6 = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab6 = scale(liab6) #6rd recurrence
liab7 = scale(g_liab) + scale(frac * rnorm(N) + (1-frac) * frail); liab7 = scale(liab7) #7th recurrence

# Binary phenotype
Y_bin = rbinom(N, 1, prob = 1/(1+exp(-liab + 1)))          # Approx 30% cases

# TTE phenotype (with 30% censoring)
Y_tte = exp(-liab); Y_tte_status = rep(1, N)
Cens = sample(Y_tte) - 0.6                                 # Approx 30% events
Y_tte_status[Y_tte > Cens] = 0; Y_tte[Y_tte > Cens] = Cens[Y_tte > Cens]
tte = cbind(1:N,Y_tte, Y_tte_status)

Y_tte_resid = coxph(Surv(Y_tte, Y_tte_status) ~ 1)$residuals

# RE phenotype (with 1 event on average) 

re = NULL

Cens = runif(N) * 2.5
for(i in 1:N){
  add = NULL
  
  add = rbind(add, c(i, exp(-liab2[i]), 1))
  add = rbind(add, c(i, exp(-liab3[i]), 1))
  add = rbind(add, c(i, exp(-liab4[i]), 1))
  add = rbind(add, c(i, exp(-liab5[i]), 1))
  add = rbind(add, c(i, exp(-liab6[i]), 1))
  add = rbind(add, c(i, exp(-liab7[i]), 1))
  
  fin = Cens[i] - cumsum(add[,2]); fin = min(Cens[i], fin[which(fin >= 0)[1]], na.rm=T)
  add = add[cumsum(add[,2]) < Cens[i],]
  
  add = rbind(add, c(i, fin, 0))
  
  re = rbind(re, add)
}

sum(re[,3] == 1)/N

Y_re_resid = rowsum(coxph(Surv(re[,2], re[,3]) ~ frailty.gaussian(re[,1]))$residuals, group = re[,1])

##################################################
####### COMPUTE ASSOCIATION STATISTICS ###########
##################################################

G_cent = apply(G, 2, scale)

# Binary
p_bin = apply(G_cent, 2, function(x) summary(glm(Y_bin ~ x))$coefficients[2,4])

# TTE
p_tte = apply(G_cent, 2, function(x) summary(lm(Y_tte_resid ~ x))$coefficients[2,4])

# RE
p_re = apply(G_cent, 2, function(x) summary(lm(Y_re_resid ~ x))$coefficients[2,4])

# Quant
p_quant = apply(G_cent, 2, function(x) summary(lm(liab ~ x))$coefficients[2,4])

ps = rbind(c('Binary',p_bin[1:22]),
           c('TTE',p_tte[1:22]),
           c('RE',p_re[1:22]),
           c('Quant',p_quant[1:22]))

all = rbind(all, ps)
print(j)

}

##########################################
#########################################
##########################################
#########################################

save = as.data.frame(all)

#save.image('Association_part.R')

#save[,-1] = apply(save[,-1],2,function(x) mean(as.numeric(x) < 0.0001))
save[,-1] = apply(save[,-1],2,function(x) qchisq(1 - as.numeric(x), df = 1))
save[save == Inf] = 70
colnames(save) = NULL

load('Prediction_part.R')

power = rbind(c('Binary',colMeans(rbind(as.matrix(save[save[,1] == 'Binary',2*(1:11)]),
                                        as.matrix(save[save[,1] == 'Binary',2*(1:11)+1])))),
              c('Time-to-event',colMeans(rbind(as.matrix(save[save[,1] == 'TTE',2*(1:11)]),
                                        as.matrix(save[save[,1] == 'TTE',2*(1:11)+1])))),
              c('Recurrent event',colMeans(rbind(as.matrix(save[save[,1] == 'RE',2*(1:11)]),
                                        as.matrix(save[save[,1] == 'RE',2*(1:11)+1])))),
              c('Quantitative',colMeans(rbind(as.matrix(save[save[,1] == 'Quant',2*(1:11)]),
                                        as.matrix(save[save[,1] == 'Quant',2*(1:11)+1])))))

colnames(power) = c('Trait',seq(from = 0, to = 1, by = 0.1))

data = pivot_longer(as.data.frame(power), -1)
data$name = as.numeric(data$name)
data$value = as.numeric(data$value)

data$Trait = as.factor(data$Trait)
data$Trait = factor(data$Trait, levels = c('Quantitative','Binary','Time-to-event','Recurrent event'))

sig = qchisq(1 - as.numeric(5e-8), df = 1)

assoc = ggplot(data, aes(name, value, col = Trait)) + geom_line() + geom_point() + xlab('SNP effect size') + ylab('Mean Statistic') + 
  geom_abline(intercept = sig, slope = 0, linetype = 'dotted') + 
  annotate("text",x = 0.3, y = sig+2, label = 'Genome-wide significance', fontface = 2) + theme_light()

####################################################################################
############################ POWER ############################
####################################################################################

# use glmnet function? For varying heritabilities?
# phenotypes liab Y_bin Y_tte_resid Y_re_resid

#install.packages('glmnet')
library(glmnet)

all = NULL

### introduction - sample sizes?
vals = 50   # will perform vals * 11 repetitions

for(j in 1:(vals*11)){
  
  N = 1000
  n_snps = 100
  her = rep(seq(from = 0, to = 1, by = 0.1), each = vals)[j]
  
  G = matrix(rbinom(N * n_snps, 2, 0.5), nrow = N)
  eff = c(rnorm(n_snps/2), rep(0, n_snps/2))
  g_liab = G %*% eff;
  
  # introduce some frailty (for RE)
  frac = 0.8
  frail = rnorm(N)
  
  liab = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab = scale(liab)    # This is already the quantitative phenotype
  
  liab2 = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab2 = scale(liab2) #2nd recurrence
  liab3 = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab3 = scale(liab3) #3rd recurrence
  liab4 = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab4 = scale(liab4) #4th recurrence
  liab5 = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab5 = scale(liab5) #5nd recurrence
  liab6 = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab6 = scale(liab6) #6rd recurrence
  liab7 = sqrt(her)*scale(g_liab) + sqrt((1-her))*scale(frac * rnorm(N) + (1-frac) * frail); liab7 = scale(liab7) #7th recurrence
  
  # Binary phenotype
  Y_bin = rbinom(N, 1, prob = 1/(1+exp(-liab + 1)))          # Approx 30% cases
  
  # TTE phenotype (with 30% censoring)
  Y_tte = exp(-liab); Y_tte_status = rep(1, N)
  Cens = sample(Y_tte) - 0.6                                 # Approx 30% events
  Y_tte_status[Y_tte > Cens] = 0; Y_tte[Y_tte > Cens] = Cens[Y_tte > Cens]
  tte = cbind(1:N,Y_tte, Y_tte_status)
  
  Y_tte_resid = coxph(Surv(Y_tte, Y_tte_status) ~ 1)$residuals
  
  # RE phenotype (with 1 event on average) 
  
  re = NULL
  
  Cens = runif(N) * 2.5
  for(i in 1:N){
    add = NULL
    
    add = rbind(add, c(i, exp(-liab2[i]), 1))
    add = rbind(add, c(i, exp(-liab3[i]), 1))
    add = rbind(add, c(i, exp(-liab4[i]), 1))
    add = rbind(add, c(i, exp(-liab5[i]), 1))
    add = rbind(add, c(i, exp(-liab6[i]), 1))
    add = rbind(add, c(i, exp(-liab7[i]), 1))
    
    fin = Cens[i] - cumsum(add[,2]); fin = min(Cens[i], fin[which(fin >= 0)[1]], na.rm=T)
    add = add[cumsum(add[,2]) < Cens[i],]
    
    add = rbind(add, c(i, fin, 0))
    
    re = rbind(re, add)
  }
  
  sum(re[,3] == 1)/N
  
  Y_re_resid = rowsum(coxph(Surv(re[,2], re[,3]) ~ 1)$residuals, group = re[,1])
  
  ##################################################
  ####### COMPUTE ASSOCIATION STATISTICS ###########
  ##################################################
  
  G_cent = apply(G, 2, scale)
  n_train = 800
  lambda = 0.5
  
  liab = scale(liab); Y_bin = scale(Y_bin); Y_tte_resid = scale(Y_tte_resid); Y_re_resid = scale(Y_re_resid)
  
  lambda = her/n_snps
  
  fit = glmnet(G_cent[1:n_train,], liab[1:n_train], family = 'gaussian', alpha = 0.5, lambda = lambda)
  out = predict(fit, newx = G_cent[-c(1:n_train),]) # make predictions
  quant = mean((liab[-c(1:n_train)] - as.numeric(out))^2)/(mean(liab[-c(1:n_train)]^2))
  
  fit = glmnet(G_cent[1:n_train,], Y_bin[1:n_train], family = 'gaussian', alpha = 0.5, lambda = lambda)
  out = predict(fit, newx = G_cent[-c(1:n_train),]) # make predictions
  bin = mean((Y_bin[-c(1:n_train)] - as.numeric(out))^2)/(mean(Y_bin[-c(1:n_train)]^2))
  
  fit = glmnet(G_cent[1:n_train,], Y_tte_resid[1:n_train], family = 'gaussian', alpha = 0.5, lambda = lambda)
  out = predict(fit, newx = G_cent[-c(1:n_train),]) # make predictions
  tte = mean((Y_tte_resid[-c(1:n_train)] - as.numeric(out))^2)/(mean(Y_tte_resid[-c(1:n_train)]^2))
  
  fit = glmnet(G_cent[1:n_train,], Y_re_resid[1:n_train], family = 'gaussian', alpha = 0.5, lambda = lambda)
  out = predict(fit, newx = G_cent[-c(1:n_train),]) # make predictions
  re = mean((Y_re_resid[-c(1:n_train)] - as.numeric(out))^2)/(mean(Y_re_resid[-c(1:n_train)]^2))
  
  all = rbind(all, c(her, quant, bin, tte, re))
  
  print(j)
}

##########################################
#########################################
##########################################
#########################################

colnames(all) = c('Her','Quantitative','Binary','Time-to-event','Recurrent event')
data = pivot_longer(as.data.frame(all), -1)
data$value = as.numeric(data$value)

comb = NULL

for(her in unique(data$Her)){
  for(name in unique(data$name)){
    row = c(her, name, mean(data$value[data$Her == her & data$name == name]))
    comb = rbind(comb, row)
  }
}

comb = as.data.frame(comb)
colnames(comb) = c('Her','name','value')
comb$Her = as.numeric(comb$Her)
comb$value = as.numeric(comb$value)

comb$name = as.factor(comb$name)
comb$name = factor(comb$name, levels = c('Quantitative','Binary','Time-to-event','Recurrent event'))

#save.image('Prediction_part.R')

pred = ggplot(comb, aes(Her, value, col = name)) + geom_line() + geom_point() + xlab('Heritability') + ylab('Mean Squared Error') + 
  geom_abline(intercept = 1, slope = 0)  +  theme_light() + labs(colour="Trait type")

###########
## Make histograms
###########

N = 1000
n_snps = 100
her = rep(seq(from = 0, to = 1, by = 0.1), each = vals)[j]

G = matrix(rbinom(N * n_snps, 2, 0.5), nrow = N)
eff = c(rnorm(n_snps/2), rep(0, n_snps/2))
g_liab = G %*% eff;

# introduce some frailty (for RE)
frac = 0.8
frail = rnorm(N)

liab = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab = scale(liab)    # This is already the quantitative phenotype
liab2 = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab2 = scale(liab2) #2nd recurrence
liab3 = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab3 = scale(liab3) #3rd recurrence
liab4 = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab4 = scale(liab4) #4th recurrence
liab5 = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab5 = scale(liab5) #5nd recurrence
liab6 = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab6 = scale(liab6) #6rd recurrence
liab7 = her*scale(g_liab) + (1-her)*scale(frac * rnorm(N) + (1-frac) * frail); liab7 = scale(liab7) #7th recurrence

# Binary phenotype
Y_bin = rbinom(N, 1, prob = 1/(1+exp(-liab + 1)))          # Approx 30% cases

# TTE phenotype (with 30% censoring)
Y_tte = exp(-liab); Y_tte_status = rep(1, N)
Cens = sample(Y_tte) - 0.6                                 # Approx 30% events
Y_tte_status[Y_tte > Cens] = 0; Y_tte[Y_tte > Cens] = Cens[Y_tte > Cens]
tte = cbind(1:N,Y_tte, Y_tte_status)

Y_tte_resid = coxph(Surv(Y_tte, Y_tte_status) ~ 1)$residuals

# RE phenotype (with 1 event on average) 

re = NULL

Cens = runif(N) * 2.5
for(i in 1:N){
  add = NULL
  
  add = rbind(add, c(i, exp(-liab2[i]), 1))
  add = rbind(add, c(i, exp(-liab3[i]), 1))
  add = rbind(add, c(i, exp(-liab4[i]), 1))
  add = rbind(add, c(i, exp(-liab5[i]), 1))
  add = rbind(add, c(i, exp(-liab6[i]), 1))
  add = rbind(add, c(i, exp(-liab7[i]), 1))
  
  fin = Cens[i] - cumsum(add[,2]); fin = min(Cens[i], fin[which(fin >= 0)[1]], na.rm=T)
  add = add[cumsum(add[,2]) < Cens[i],]
  
  add = rbind(add, c(i, fin, 0))
  
  re = rbind(re, add)
}

sum(re[,3] == 1)/N

Y_re_resid = rowsum(coxph(Surv(re[,2], re[,3]) ~ 1)$residuals, group = re[,1])

### make hists

pheno = cbind(liab, Y_bin, Y_tte_resid, Y_re_resid)
pheno = as.data.frame(pheno)
colnames(pheno) = c('Quantitative','Binary',"Time","Recurrent")

plots = list()

for(i in 1:4){
  plots[[i]] = ggplot(pheno, aes_string(x = colnames(pheno)[i])) + geom_histogram() + xlab('Phenotype') + ylab('Frequency') + 
  geom_abline(intercept = 1, slope = 0)  +  theme_classic() + ggtitle(colnames(pheno)[i]) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10)
  )
}

plots[[1]] = plots[[1]] + annotate('text', x = -2.5,y = 45, label  = 'A', size = 10) + xlim(c(-3.5,3.5))
plots[[2]] = plots[[2]] + scale_x_continuous(breaks = c(0,1))
plots[[3]] = plots[[3]] + ggtitle('Time-to-event')
plots[[4]] = plots[[4]] + ggtitle('Recurrent event')

library(gridExtra)
library(grid)

top <- textGrob(expression("Phenotype distribution"), gp = gpar(fontsize = 17), hjust = 0.45)
hist = grid.arrange(plots[[1]], plots[[2]],plots[[3]],plots[[4]], nrow = 2,
                    top = top)

save.image('Prediction_part.R')

############## Combine all plots!!! ###############

setwd("//umcsanfsclp01/hev_genepi/genepi/Jasper/Papers0_Thesis/Scripts")

load('Prediction_part.R')

library(ggplot2)
library(gridExtra)
library(grid)

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

pred = pred + theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16)
) + annotate('text', x = 0.1, y = 0.9, label = 'C', size = 10)
assoc = assoc + theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
) + annotate('text', x = 0.1, y = 50, label = 'B', size = 10)

legend = get_only_legend(pred)




assoc = assoc + theme(legend.position = 'none',
                      plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))
pred = pred + theme(legend.position = 'none',
                    plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm")) + xlab('Heritability of Linear Predictor')

top <- textGrob(expression("Statistical power"), gp = gpar(fontsize = 17), hjust = 0.4)
assoc_new = grid.arrange(assoc, nrow = 1,
                    top = top)

top <- textGrob(expression("PRS accuracy"), gp = gpar(fontsize = 17), hjust = 0.4)
pred_new = grid.arrange(pred, nrow = 1,
                         top = top)

plots = list(hist, assoc_new, pred_new, legend)

png("//umcsanfsclp01/hev_genepi/genepi/Jasper/Papers0_Thesis/Figures/trait_comparison.png",
    width = 15, height = 5, units = 'in', res = 1200)

grid.arrange(grobs = plots, nrow = 1, widths = c(0.3,0.3,0.3,0.15))

dev.off()



