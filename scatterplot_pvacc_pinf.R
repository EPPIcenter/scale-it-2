library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(simstudy)
library(readr)
library(lubridate)
library(tidyverse)
library(patchwork)

#############################
## Figure: p_vacc vs p_inf ##
#############################

load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v2.RData")
load("~/../Box/PHI_scaleit/SCALEIT_2.0/labels.RData")

## Look at 2D
summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "p_inf") | str_starts(parameter, "p_vacc\\[")) %>%
  separate(col=parameter, into=c("tmp","parameter","demog_group")) %>%
  select(-tmp) %>%
  mutate(parameter=paste0("p_", parameter)) %>%
  mutate(demog_group = as.numeric(demog_group)) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  ## Update labels
  ungroup() %>%
  mutate(race = ifelse(race==" Asian", "Asian", race)) %>%
  mutate(race = ifelse(race==" White", "White", race)) %>%
  mutate(race = ifelse(race==" Black or African American", "Black", race)) %>%
  mutate(race = ifelse(race==" Hispanic or Latino", "Latinx", race)) %>%
  as.data.frame() %>%
  ## Rename
  rename(Age = age) %>%
  rename(`Race/Ethnicity` = race) %>%
  pivot_wider(names_from=parameter, values_from=mean:Rhat) %>%
  ggplot() +
  geom_abline(linetype="dotted") +
  geom_point(aes(x=mean_p_inf, y=mean_p_vacc, colour=Age, shape=`Race/Ethnicity`), size=4) +
  geom_linerange(aes(x=mean_p_inf, y=mean_p_vacc, xmin=`2.5%_p_inf`, xmax=`97.5%_p_inf`, colour=Age), size=1) +
  geom_linerange(aes(x=mean_p_inf, y=mean_p_vacc, ymin=`2.5%_p_vacc`, ymax=`97.5%_p_vacc`, colour=Age), size=1) +
  theme_bw(base_size=14) +
  scale_shape_manual(values=c(15, 16, 17, 21)) +
  xlab("Probability of infection") + ylab("Probability of vaccination") +
  xlim(0,0.45) + ylim(0,0.45)

#################
## Figure: RRR ##
#################

load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v2.RData")
load("~/../Box/PHI_scaleit/SCALEIT_2.0/labels.RData")

## Pull out values to tabulate
summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "RRR_vs_white")) %>%
  separate(col=parameter, into=c("tmp1","tmp2","tmp3","demog_group")) %>%
  select(-tmp1, -tmp2, -tmp3) %>%
  mutate(parameter="RRR_vs_white_of_that_age") %>%
  mutate(demog_group = as.numeric(demog_group)) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  as.data.frame() %>%
  mutate(race = ifelse(race==" Asian", "Asian", race)) %>%
  mutate(race = ifelse(race==" White", "White", race)) %>%
  mutate(race = ifelse(race==" Black or African American", "Black", race)) %>%
  mutate(race = ifelse(race==" Hispanic or Latino", "Latinx", race)) %>%
  mutate(age = ifelse(age=="18-64 ", "18-64", "65+")) %>%
  mutate(label_for_plot2 = paste0(age, " & ", race, " (n=", n, ")")) %>%
  ggplot() +
  geom_vline(xintercept=1, colour=2) +
  facet_wrap(.~age, ncol=1, scales="free_y") +
  geom_pointrange(aes(y=label_for_plot2, x=`50%`, xmin=`2.5%`, xmax=`97.5%`)) +
  theme_bw() +
  xlab("Relative risk ratio") +
  ylab("") +
  scale_x_log10(breaks=scales::breaks_log(6)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())
