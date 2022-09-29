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

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

a <- read_csv(here("data","cumulative_risk.csv"))

a %>%
  ggplot(.) +
  scale_y_continuous(expand = c(0,0), limits = c(0, .5), breaks = seq(0,.5,.1)) + 
  scale_x_continuous(expand = c(0,0), limits = c(0,90)) + 
  ylab("Conceptions Subset") + 
  xlab("Week on Study") + 
  geom_step(aes(x = time, 
                y = cum_risk, 
                color = factor(exposure)),
            alpha=.25) +
  facet_wrap(~method)

a0 <- a %>% 
        group_by(index, exposure, method) %>% 
        filter(exposure==0, 
               row_number()==n()) %>% 
  mutate(cum_risk0 = cum_risk) %>% 
  ungroup() %>% 
  select(cum_risk0, method, index, -exposure)

a1 <- a %>% 
  group_by(index, exposure, method) %>% 
  filter(exposure==1, 
         row_number()==n())  %>% 
  mutate(cum_risk1 = cum_risk) %>% 
  ungroup() %>% 
  select(cum_risk1, method, index, -exposure) 

table(a0$method)
table(a1$method)

a_rd <- left_join(a0,a1) %>% mutate(risk_difference = cum_risk1 - cum_risk0)

a_rd %>% ggplot(.) +
  geom_histogram(aes(risk_difference)) +
  facet_wrap(~method)

a_rd %>% group_by(method) %>% summarize(mean_res = mean(risk_difference),
                                        median_res = median(risk_difference),
                                        se_res = sd(risk_difference))