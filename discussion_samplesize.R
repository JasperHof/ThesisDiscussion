
library(readxl)
library(survminer)
library(survival)

####################
### introduction - sample sizes?
###################

N = 4741
all_p = NULL


eff = 0.19 # recurrence
eff = 0.32 # recurrence

2066/N  # 43% recurrence rate
696/N   # 15% recurrence rate

for(i in 1:100){

g = rbinom(N, 2, 0.30)
g_liab = g * eff;
liab = g_liab + scale(rnorm(N)); liab = scale(liab)    # This is already the quantitative phenotype

# TTE phenotype (with 30% censoring)
Y_tte = exp(-liab); Y_tte_status = rep(1, N)

#Cens = sample(Y_tte) - 0.2                              # Approx 43% events
Cens = sample(Y_tte) - 1.7                            # Approx 43% events

Y_tte_status[Y_tte > Cens] = 0; Y_tte[Y_tte > Cens] = Cens[Y_tte > Cens]
tte = cbind(1:N,Y_tte, Y_tte_status)

all_p = c(all_p, summary(coxph(Surv(Y_tte, Y_tte_status) ~ g))$coefficients[5])

print(i)

}

mean(all_p < 5e-8)

exp(0.19)
exp(0.32)
