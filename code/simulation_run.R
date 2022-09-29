#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## NOTE: We can run this program locally from the command line, and thus leave the package 
## installation as is. Once we run on cluster, need to change library location.

## where is the library?
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","skimr","here","survival","VGAM","cmprsk")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

# baseline hazard: Weibull

# N = sample size    
# lmbd = scale parameter for weibull
# ro = shape parameter for weibull
# bta = regression parameter
# rteC = rate parameter for censoring 
# admnC = administrative censoring time

seed = as.numeric(args[1])
N = as.numeric(args[2])
lmbd <- c(.005, .005, .005, .00525)
ro <- c(.75, .65, .85, .6)
bta <- c(1.2,1.7,.7,.5)
rte <- c(.05, .04, .06, .025)
admnC <- 90

weibullSim <- function(N, lambda, rho, beta, rateC, adminC, seed){
  
  set.seed(seed)
  
  x <- rbinom(n = N, 1, .5)
  # inverse transform sampling, modified from Bender et al 2005, p 1717, Table II
  ## uniform distribution: one for each competing event
  u <- cbind(runif(n = N),runif(n = N),runif(n = N),runif(n = N))
  ## four competing time to events
  T_outcome1 <- (- log(u[,1]) / (lambda[1] * exp(x * beta[1])))^(1 / rho[1])
  T_outcome2 <- (- log(u[,2]) / (lambda[2] * exp(x * beta[2])))^(1 / rho[2])
  T_outcome3 <- (- log(u[,3]) / (lambda[3] * exp(x * beta[3])))^(1 / rho[3])
  T_outcome4 <- (- log(u[,4]) / (lambda[4] * exp(x * beta[4])))^(1 / rho[4])
  ## exponential censoring
  C <- rexp(n=N, rate=rateC)
  ## find the min of the time to events to determine which is "seen" first
  T_outcome <- pmin(pmin(pmin(T_outcome1,T_outcome2),T_outcome3),T_outcome4)
  ## does censoring occur before the outcome?
  time <- pmin(T_outcome, C)
  ## code the event indicator, 0 for censoring
  event <-  1*(T_outcome1==T_outcome)*as.numeric(T_outcome <= C) +
             2*(T_outcome2==T_outcome)*as.numeric(T_outcome <= C) +
              3*(T_outcome3==T_outcome)*as.numeric(T_outcome <= C) +
               4*(T_outcome4==T_outcome)*as.numeric(T_outcome <= C)
  ## define the end of study follow up
  time <- ifelse(time>adminC, adminC, time)
  event <- ifelse(time>adminC, 0, event)
  # data set

  sim_dat <- data.frame(id=1:N,
                        time=time,
                        event=event,
                        exposure=x)

  ## Aalen Johansen via survfit()
  aj_mod1 <- survfit(
    Surv(time = time, 
         event = factor(event)) ~ 1, 
    data = subset(sim_dat, exposure==1))
  aj_mod0 <- survfit(
    Surv(time = time, 
         event = factor(event)) ~ 1, 
    data = subset(sim_dat, exposure==0))
  
  aj_dat1 <- tibble(time = aj_mod1$time, 
                    cum_risk = aj_mod1$pstate[, 2],
                    exposure=1, 
                    method="Aalen-Johansen",
                    index = seed)
  aj_dat0 <- tibble(time = aj_mod0$time, 
                    cum_risk = aj_mod0$pstate[,2],
                    exposure=0, 
                    method="Aalen-Johansen",
                    index = seed)
  
  cum_risk_AJ <- rbind(aj_dat1, aj_dat0)
  
  # multinomial GLM via VGAM
  vglm_func <- function(data, index){
    a <- data[index,] %>% 
      mutate(ID = 1:n(),
             stop_rounded = ceiling(time)) %>% 
      uncount(stop_rounded) %>% 
      group_by(ID) %>% 
      mutate(counter = 1,
             time_var = cumsum(counter),
             last_id = !duplicated(ID, fromLast = T),
             outcome = event*last_id) %>% 
      ungroup(ID) %>% 
      select(ID, time_var, outcome, time, exposure, event)
    
    # model time-to-events via multinomial
    mod1 <- vglm(outcome ~ ns(time_var, knots=c(10,15,20,25,30,40)),  
                 data=subset(a, exposure==1), 
                 family=multinomial(refLevel = 1))
    mod0 <- vglm(outcome ~ ns(time_var, knots=c(10,15,20,25,30,40)),  
                 data=subset(a, exposure==0), 
                 family=multinomial(refLevel = 1))
    
    # predict risks for outcome of interest
    mu_11 <- tibble(a, mu_11 = predict(mod1, 
                                       newdata = a, 
                                       type = "response")[, 2])
    mu_10 <- tibble(a, mu_10 = predict(mod0, 
                                       newdata = a, 
                                       type = "response")[, 2])
    
    # average risks across sample at each time point
    mu_11 <- mu_11 %>% 
      group_by(time_var) %>% 
      summarize(mean_mu_11 = mean(mu_11)) 
    mu_10 <- mu_10 %>% 
      group_by(time_var) %>% 
      summarize(mean_mu_10 = mean(mu_10))
    
    # cumulate risks over time
    mu_11 <- mu_11 %>% 
      mutate(cum_risk = cumsum(mean_mu_11))
    mu_10 <- mu_10 %>% 
      mutate(cum_risk = cumsum(mean_mu_10))
  
    results <- list(data.frame(time = mu_10$time_var,
                               cum_risk0 = mu_10$cum_risk),
                    data.frame(time = mu_11$time_var,
                               cum_risk1 = mu_11$cum_risk))
    
    return(results)
  }
  
  res <- vglm_func(sim_dat)
  
  res0 <- res[[1]] %>% mutate(cum_risk=cum_risk0, exposure=0, method = "VGLM", index = seed) %>% select(time, cum_risk, exposure, method, index)
  res1 <- res[[2]] %>% mutate(cum_risk=cum_risk1, exposure=1, method = "VGLM", index = seed) %>% select(time, cum_risk, exposure, method, index)
  
  cum_risk_VGLM <- rbind(res0,res1)
  
  # Gray's subdistribution function estimator
  res_gray <- cuminc(ftime=sim_dat$time, 
                     fstatus=sim_dat$event,
                     group=sim_dat$exposure,
                     strata=sim_dat$exposure)
  
  res_dat0 <- data.frame(time = res_gray$`0 2`$time,
                         cum_risk = res_gray$`0 2`$est) %>% 
    mutate(exposure = 0, method = "Gray", index=seed)
  res_dat1 <- data.frame(time = res_gray$`1 2`$time,
                         cum_risk = res_gray$`1 2`$est) %>% 
    mutate(exposure = 1, method = "Gray", index=seed)
  
  cum_risk_Gray <- rbind(res_dat0, res_dat1)
  
  cumulative_risk <- rbind(cum_risk_AJ,
                           cum_risk_VGLM,
                           cum_risk_Gray)
  
  return(cumulative_risk)
  
}

sim_res <- lapply(1:seed, 
                  function(x)
                  weibullSim(N = N, 
                             lambda = lmbd, 
                             rho = ro, 
                             beta = bta, 
                             rateC = rte, 
                             adminC = admnC,
                             seed = x)
)

write_csv(here("data","cumulative_risk.csv"))