
library(readxl)
library(survminer)
library(survival)

####################
### introduction - sample sizes?
###################

N = 3400
all_p = NULL


#eff = 0.21 # recurrence
#eff = 0.35 # progression

#eff = log(1.3) # recurrence (old)
eff = 0.43 # progression (old)

2145/N  # 43% recurrence rate
#696/N   # 15% recurrence rate

for(i in 1:100){

g = rbinom(N, 2, 0.20)
g_liab = g * eff;
liab = g_liab + scale(rnorm(N)); liab = scale(liab)    # This is already the quantitative phenotype

# TTE phenotype (with 30% censoring)
Y_tte = exp(-liab); Y_tte_status = rep(1, N)

#Cens = sample(Y_tte) - 0.2                              # Approx 43% events
Cens = sample(Y_tte) - 1.7                            # Approx 15% events

# old
#Cens = sample(Y_tte) - 0.1                              # Approx 43% events

Y_tte_status[Y_tte > Cens] = 0; Y_tte[Y_tte > Cens] = Cens[Y_tte > Cens]
tte = cbind(1:N,Y_tte, Y_tte_status)

all_p = c(all_p, summary(coxph(Surv(Y_tte, Y_tte_status) ~ g))$coefficients[5])

print(i)

}

mean(all_p < 5e-8)

exp(eff)


exp(0.19)
exp(0.32)

exp(0.43)

### multiple recurrences

library(simrec)
library(survival)

N = 5009
all_p = NULL
eff = 0.18 # recurrence

4237/N  # 85%% recurrence rate

dist.x <- c("binomial", "binomial")
par.x <- list(0.2, 0.2)
par.rec <- c(1,2)
par.eff = c(eff, eff)

for(i in 1:50){
  
  g = rbinom(N, 2, 0.30)
  g_liab = g * eff;
  liab = g_liab + scale(rnorm(N)); liab = scale(liab)    # This is already the quantitative phenotype
  
  
  data = simrec(N, 0, 1.54, dist.x = dist.x, par.x = par.x,
         beta.x = par.eff, dist.z = "gamma", par.z = 0.2, dist.rec = "weibull", par.rec = par.rec, pfree = 0,
         dfree = 0)
  
  data$G = data$x.V1 + data$x.V2
  fit = coxph(Surv(start, stop, status) ~ G + frailty.gamma(id), data = data)
  
  p = summary(fit)$coefficients[1,6]
  
  all_p = c(all_p, p)
  
  print(mean(all_p < 5e-8))
  print(i)
}


exp(0.20)

mean(all_p < 5e-8)

exp(0.19)
exp(0.32)


###############

# sample size calculations nbcs-2 urolife
library(survival)

N = 1630
all_p = NULL

#eff = 0.16 # recurrence
eff = 0.27 # progression


for(i in 1:100){
  
  g = rbinom(N, 2, 0.20)
  g_liab = g * eff;
  liab = g_liab + scale(rnorm(N)); liab = scale(liab)    # This is already the quantitative phenotype
  
  # TTE phenotype (with 30% censoring)
  Y_tte = exp(-liab); Y_tte_status = rep(1, N)
  
  #Cens = sample(Y_tte) - 0.2                              # Approx 43% events
  Cens = sample(Y_tte) - 1.7                            # Approx 15% events
  
  # old
  #Cens = sample(Y_tte) - 0.1                              # Approx 43% events
  
  Y_tte_status[Y_tte > Cens] = 0; Y_tte[Y_tte > Cens] = Cens[Y_tte > Cens]
  tte = cbind(1:N,Y_tte, Y_tte_status)
  
  all_p = c(all_p, summary(coxph(Surv(Y_tte, Y_tte_status) ~ g))$coefficients[5])
  
  print(i)
  
}

mean(all_p < 5e-2)

exp(0.16) # recurrence
exp(0.27) # progression
