library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(simstudy)
library(readr)
library(lubridate)
library(tidyverse)
library(patchwork)

source("~/../Box/PHI_scaleit/SCALEIT_2.0/scripts/Stan_functions_v2.R")

## Read in data
dat_SCALEIT <- read_csv("~/../Box/PHI_scaleit/SCALEIT_2.0/reshaped_redux.csv")

## Look at the raw univariate results
dat_SCALEIT %>%
  group_by(seropos_Roche) %>%
  tally()

dat_SCALEIT %>%
  group_by(seropos_Vitros) %>%
  tally()

dat_SCALEIT %>%
  group_by(seropos_Vitros,seropos_Roche) %>%
  tally()

dat_SCALEIT %>%
  group_by(age_groups) %>%
  tally()

dat_SCALEIT %>%
  group_by(DateCollected) %>%
  tally()

## Plot
dat_SCALEIT %>%
  ## Remove NAs
  filter(!is.na(seropos_Vitros) & !is.na(seropos_Roche)) %>%
  filter(!is.na(age_groups)) %>%
  mutate(age_groups = factor(age_groups)) %>%
  ggplot(aes(x=SC_Roche, y=SC_Vitros)) +
  geom_point(size=3, alpha=0.3, colour="cornflowerblue") +
  theme_bw() +
  ggtitle("SCALE-IT 2.0 (February 2021)") +
  scale_x_continuous(trans=scales::pseudo_log_trans(base=10), breaks=c(0,10^c(0:4)), limits=c(0,10^4)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10), breaks=c(0,10^c(0:4)), limits=c(0,10^4)) +
  facet_wrap(.~age_groups, labeller="label_both") +
  theme(
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank())

## Calculate the weighted mean date
## Use this to decide which date to use for reported vacc coverage
dat_SCALEIT %>%
  group_by(DateCollected) %>%
  tally() %>%
  ungroup() %>%
  summarise(median = weighted.mean(DateCollected, w=n))

## Separate the raw univariate data by age (2 groups) and race/ethnicity (4 groups)
dat_SCALEIT %>%
  ## Remove NAs
  filter(!is.na(seropos_Vitros) & !is.na(seropos_Roche)) %>%
  filter(!is.na(age_groups)) %>%
  ## Remove youngest age group, since ineligible for vaccination
  filter(age_groups!="0-17") %>%
  mutate(age_groups_binary = ifelse(age_groups=="65+", "65+", "18-64")) %>%
  ## Remove "Other" race
  filter(PtRace!="Other") %>%
  mutate(age_race = paste0(age_groups_binary, " & ", PtRace)) %>%
  group_by(age_race, seropos_Vitros, seropos_Roche) %>%
  tally() %>%
  ungroup() %>%
  ## Complete missing combos for Vitros & Roche by age_race
  complete(seropos_Vitros, seropos_Roche, nesting(age_race), fill=list(n=0)) %>%
  arrange(age_race, desc(seropos_Vitros), desc(seropos_Roche)) %>%
  ## Remove Roche positive, Vitros negative results
  filter((seropos_Vitros==0 & seropos_Roche==1)==FALSE) %>%
  select(age_race, seropos_Vitros, seropos_Roche, n) %>%
  as.data.frame() -> dat_for_Stan

## -----------------------------------------------------------

######################
## Implement priors ##
######################

## Data from: https://data.sfgov.org/COVID-19/COVID-19-Vaccine-Doses-Given-to-San-Franciscans-by/xjh5-h442
dat_vacc_raw <- read_csv("~/../Box/PHI_scaleit/SCALEIT_2.0/COVID-19_Vaccine_Doses_Given_to_San_Franciscans_by_Demographics_Over_Time_2021-06-17.csv")

## Rename demographic groups
dat_vacc_raw %>%
  mutate(DEMOGRAPHIC_SUBGROUP = ifelse(DEMOGRAPHIC_SUBGROUP=="Hispanic or Latino/a, all races", "Hispanic or Latino", DEMOGRAPHIC_SUBGROUP)) -> dat_vacc_raw

## Subset & pull out the demographic data
dat_vacc_raw %>%
  filter(DEMOGRAPHIC_SUBGROUP %in% c("Asian", "Black or African American", "Hispanic or Latino", "White")) %>%
  group_by(AGE_GROUP, DEMOGRAPHIC_SUBGROUP) %>%
  summarise(n=mean(SUBGROUP_POPULATION)) %>%
  ungroup() -> dat_demog

dat_demog %>%
  add_row(AGE_GROUP="12-64", DEMOGRAPHIC_SUBGROUP="Asian", n=dat_demog$n[1]-dat_demog$n[5]) %>%
  add_row(AGE_GROUP="12-64", DEMOGRAPHIC_SUBGROUP="Black or African American", n=dat_demog$n[2]-dat_demog$n[6]) %>%
  add_row(AGE_GROUP="12-64", DEMOGRAPHIC_SUBGROUP="Hispanic or Latino", n=dat_demog$n[3]-dat_demog$n[7]) %>%
  add_row(AGE_GROUP="12-64", DEMOGRAPHIC_SUBGROUP="White", n=dat_demog$n[4]-dat_demog$n[8]) -> dat_demog

dat_demog %>%
  filter(AGE_GROUP=="12-64") -> dat_demog_12_64

## Calculate doses among 65+ year olds by date and race/ethnicity
dat_vacc_raw %>%
  filter(ADMINISTERING_PROVIDER_TYPE=="All Providers") %>%
  filter(AGE_GROUP=="65+") %>%
  filter(DEMOGRAPHIC_GROUP=="Race/Ethnicity") %>%
  mutate(date_admin = mdy(word(DATE_ADMINISTERED, 1))) %>%
  select(date_admin, AGE_GROUP, DEMOGRAPHIC_SUBGROUP, DEMOGRAPHIC_SUBGROUP_SORT_ORDER, NEW_RECIPIENTS, CUMULATIVE_RECIPIENTS, SUBGROUP_POPULATION, AGE_GROUP_POPULATION) %>%
  arrange(date_admin, DEMOGRAPHIC_SUBGROUP_SORT_ORDER) %>%
  ## Limit to certain groups
  filter(DEMOGRAPHIC_SUBGROUP %in% c("Asian", "Black or African American", "Hispanic or Latino", "White")) -> dat_vacc_65plus_raceEth

## Calculate doses among 12-64 year olds by date and race/ethnicity
dat_vacc_raw %>%
  filter(ADMINISTERING_PROVIDER_TYPE=="All Providers") %>%
  filter(AGE_GROUP=="12+") %>%
  filter(DEMOGRAPHIC_GROUP=="Race/Ethnicity") %>%
  mutate(date_admin = mdy(word(DATE_ADMINISTERED, 1))) %>%
  select(date_admin, AGE_GROUP, DEMOGRAPHIC_SUBGROUP, DEMOGRAPHIC_SUBGROUP_SORT_ORDER, NEW_RECIPIENTS, CUMULATIVE_RECIPIENTS, SUBGROUP_POPULATION, AGE_GROUP_POPULATION) %>%
  arrange(date_admin, DEMOGRAPHIC_SUBGROUP_SORT_ORDER) %>%
  ## Limit to certain groups
  filter(DEMOGRAPHIC_SUBGROUP %in% c("Asian", "Black or African American", "Hispanic or Latino", "White")) %>%
  rename(NEW_RECIPIENTS_12plus = NEW_RECIPIENTS) %>%
  rename(CUMULATIVE_RECIPIENTS_12plus = CUMULATIVE_RECIPIENTS) %>%
  select(-DEMOGRAPHIC_SUBGROUP_SORT_ORDER, -SUBGROUP_POPULATION, -AGE_GROUP_POPULATION) %>%
  left_join(dat_vacc_65plus_raceEth %>%
              rename(NEW_RECIPIENTS_65plus = NEW_RECIPIENTS) %>%
              rename(CUMULATIVE_RECIPIENTS_65plus = CUMULATIVE_RECIPIENTS) %>%
              select(-AGE_GROUP, -DEMOGRAPHIC_SUBGROUP_SORT_ORDER, -SUBGROUP_POPULATION, -AGE_GROUP_POPULATION), by=c("date_admin", "DEMOGRAPHIC_SUBGROUP")) %>%
  mutate(NEW_RECIPIENTS_65plus = ifelse(is.na(NEW_RECIPIENTS_65plus), 0, NEW_RECIPIENTS_65plus)) %>%
  mutate(CUMULATIVE_RECIPIENTS_65plus = ifelse(is.na(CUMULATIVE_RECIPIENTS_65plus), 0, CUMULATIVE_RECIPIENTS_65plus)) %>%
  mutate(NEW_RECIPIENTS_12to64 = NEW_RECIPIENTS_12plus - NEW_RECIPIENTS_65plus) %>%
  mutate(CUMULATIVE_RECIPIENTS_12to64 = CUMULATIVE_RECIPIENTS_12plus - CUMULATIVE_RECIPIENTS_65plus) %>%
  mutate(AGE_GROUP="12-64") %>%
  select(-NEW_RECIPIENTS_12plus, -CUMULATIVE_RECIPIENTS_12plus, -NEW_RECIPIENTS_65plus, -CUMULATIVE_RECIPIENTS_65plus) %>%
  rename(NEW_RECIPIENTS = NEW_RECIPIENTS_12to64) %>%
  rename(CUMULATIVE_RECIPIENTS = CUMULATIVE_RECIPIENTS_12to64) %>%
  mutate(SUBGROUP_POPULATION = NA) %>%
  mutate(SUBGROUP_POPULATION = ifelse(DEMOGRAPHIC_SUBGROUP=="Asian", dat_demog_12_64$n[1], SUBGROUP_POPULATION)) %>%
  mutate(SUBGROUP_POPULATION = ifelse(DEMOGRAPHIC_SUBGROUP=="Black or African American", dat_demog_12_64$n[2], SUBGROUP_POPULATION)) %>%
  mutate(SUBGROUP_POPULATION = ifelse(DEMOGRAPHIC_SUBGROUP=="Hispanic or Latino", dat_demog_12_64$n[3], SUBGROUP_POPULATION)) %>%
  mutate(SUBGROUP_POPULATION = ifelse(DEMOGRAPHIC_SUBGROUP=="White", dat_demog_12_64$n[4], SUBGROUP_POPULATION)) -> dat_vacc_12_64_raceEth

## -----------------------------------------------------------

########################################################################
## Join with SCALE-IT denominators to generate prior for vaccinations ##
########################################################################

## SCALE-IT pop by age & race
dat_SCALEIT %>%
  ## Remove NAs
  filter(!is.na(seropos_Vitros) & !is.na(seropos_Roche)) %>%
  filter(!is.na(age_groups)) %>%
  ## Remove youngest age group, since ineligible for vaccination
  filter(age_groups!="0-17") %>%
  ## Remove "Other" race
  filter(PtRace!="Other") %>%
  mutate(age_groups_binary = ifelse(age_groups=="65+", "65+", "18-64")) %>%
  mutate(age_race = paste0(age_groups_binary, " & ", PtRace)) %>%
  ## Remove Roche positive, Vitros negative results
  filter((seropos_Vitros==0 & seropos_Roche==1)==FALSE) %>%
  group_by(age_groups_binary, PtRace, age_race) %>%
  tally() %>%
  ungroup() %>%
  rename(n_in_SCALEIT = n) -> to_join_age_race

## Generate empirical distributions for beta priors
dat_vacc_12_64_raceEth %>%
  filter(date_admin == ymd("2021-02-10")-weeks(3)) %>% ## Weighted mean date
  group_by(date_admin, DEMOGRAPHIC_SUBGROUP) %>%
  summarise(
    n_cum_recipients = sum(CUMULATIVE_RECIPIENTS),
    n_denom = sum(SUBGROUP_POPULATION),
    prop_vacc_emp = n_cum_recipients / n_denom) %>%
  ungroup() %>%
  left_join(to_join_age_race %>% filter(age_groups_binary=="18-64"), by=c("DEMOGRAPHIC_SUBGROUP"="PtRace")) %>%
  filter(!is.na(n_in_SCALEIT)) -> for_sim_prior_12_64_race

dat_vacc_65plus_raceEth %>%
  filter(DEMOGRAPHIC_SUBGROUP %in% c("Under 12", "12-17") == FALSE) %>%
  filter(date_admin == ymd("2021-02-10")-weeks(3)) %>% ## Weighted mean date
  group_by(date_admin, DEMOGRAPHIC_SUBGROUP) %>%
  summarise(
    n_cum_recipients = sum(CUMULATIVE_RECIPIENTS),
    n_denom = sum(SUBGROUP_POPULATION),
    prop_vacc_emp = n_cum_recipients / n_denom) %>%
  ungroup() %>%
  left_join(to_join_age_race %>% filter(age_groups_binary=="65+"), by=c("DEMOGRAPHIC_SUBGROUP"="PtRace")) %>%
  filter(!is.na(n_in_SCALEIT)) -> for_sim_prior_65plus_race

for_sim_prior_age_race <- bind_rows(for_sim_prior_12_64_race, for_sim_prior_65plus_race)

## Generate observations from hypergeometric distributions for priors
set.seed(234324)
n_draws <- 10000

y_obs <- matrix(NA, nrow(for_sim_prior_age_race), n_draws)

for(i in 1:nrow(for_sim_prior_age_race)) {
  
  y_obs[i,] <- rhyper(n_draws, m=for_sim_prior_age_race$n_cum_recipients[i], n=for_sim_prior_age_race$n_denom[i]-for_sim_prior_age_race$n_cum_recipients[i], k=for_sim_prior_age_race$n_in_SCALEIT[i])/for_sim_prior_age_race$n_in_SCALEIT[i]

}

## No zeros
y_obs[which(y_obs==0)] <- 1e-6

## -----------------------------------------------------------

###########################
## Implement Stan models ##
###########################

## Fit the model
fit_stan <- stan(
  model_code=model_stan_Ngroup_2pvacc_priorData,
  data=list(
    n_groups=length(unique(dat_for_Stan$age_race)),
    A_plus_C=dat_for_Stan %>% filter(seropos_Vitros==1 & seropos_Roche==1) %>% select(n) %>% pull(),
    B=dat_for_Stan %>% filter(seropos_Vitros==1 & seropos_Roche==0) %>% select(n) %>% pull(),
    D=dat_for_Stan %>% filter(seropos_Vitros==0 & seropos_Roche==0) %>% select(n) %>% pull(),
    n_draws=n_draws,
    y_obs=y_obs))

#save(dat_for_Stan, fit_stan, file="~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v1.RData")
save(dat_for_Stan, fit_stan, file="~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v2.RData")

## Fit the model again, now making White the reference group
fit_stan <- stan(
  model_code=model_stan_Ngroup_2pvacc_priorData_with_RRR_vs_white,
  data=list(
    n_groups=length(unique(dat_for_Stan$age_race)),
    A_plus_C=dat_for_Stan %>% filter(seropos_Vitros==1 & seropos_Roche==1) %>% select(n) %>% pull(),
    B=dat_for_Stan %>% filter(seropos_Vitros==1 & seropos_Roche==0) %>% select(n) %>% pull(),
    D=dat_for_Stan %>% filter(seropos_Vitros==0 & seropos_Roche==0) %>% select(n) %>% pull(),
    n_draws=n_draws,
    y_obs=y_obs))

#save(dat_for_Stan, fit_stan, file="~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v1.RData")
#save(dat_for_Stan, fit_stan, file="~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v2.RData")
save(dat_for_Stan, fit_stan, file="~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v3.RData")

## -----------------------------------------------------------

################################
## Look at parameters and RRs ##
################################

## Load in model fit
#load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v1.RData")
load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v2.RData")

## Get labels for plotting
dat_for_Stan %>%
  group_by(age_race) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  mutate(label_for_plot = paste0(age_race, " (n=", n, ")")) %>%
  as.data.frame() -> dat_label_age_race

## Summary
color_scheme_set("brightblue")
fit_stan

mcmc_trace(fit_stan) + theme_bw()

mcmc_hist(fit_stan, pars=vars(`p_vacc_prior[1]`:`p_vacc_prior[8]`)) +
  theme_bw()
  
mcmc_intervals(fit_stan, pars=vars(`mu_p_vacc[1]`:`mu_p_vacc[8]`))

## Get labels for plotting
dat_for_Stan %>%
  group_by(age_race) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  mutate(label_for_plot = paste0(age_race, " (n=", n, ")")) %>%
  as.data.frame() -> dat_label_age_race

mcmc_intervals(fit_stan, pars=vars(`p_vacc[1]`:`p_vacc[8]`, `p_inf[1]`:`p_inf[8]`)) +
  ggtitle("Probability of vaccination vs. infection (model with separate p_vacc by infection status & priors)") +
  xlim(0,1)

mcmc_intervals(fit_stan, pars=vars(contains("ratio_p_vacc_p_inf")), prob_outer=0.95) +
  scale_x_log10(breaks=scales::breaks_log(6)) +
  scale_y_discrete(labels=dat_label_age_race$label_for_plot) +
  ggtitle("Relative risk of vaccination vs. infection (model with separate p_vacc by infection status & priors)") +
  theme_bw(base_size=12) +
  geom_vline(xintercept=1, colour="red")

## New labels
## Labels - age (either fine or coarse) & race
dat_for_Stan %>%
  group_by(age_race) %>%
  tally() %>%
  ungroup() %>%
  mutate(label_for_plot_inf = paste0(age_race, "; infected")) %>%
  select(label_for_plot_inf) %>%
  pull() -> label_inf

dat_for_Stan %>%
  group_by(age_race) %>%
  tally() %>%
  ungroup() %>%
  mutate(label_for_plot_uninf = paste0(age_race, "; uninfected")) %>%
  select(label_for_plot_uninf) %>%
  pull() -> label_uninf

mcmc_intervals(fit_stan, pars=vars(`p_vacc_given_inf[1]`:`p_vacc_given_inf[8]`, `p_vacc_given_uninf[1]`:`p_vacc_given_uninf[8]`), prob_outer=0.95) +
  scale_y_discrete(labels=c(label_inf, label_uninf)) +
  ggtitle("p_vacc") +
  theme_bw(base_size=12)

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
  as.data.frame() %>%
  pivot_wider(names_from=parameter, values_from=mean:Rhat) %>%
  ggplot() +
  geom_abline(linetype="dotted") +
  geom_point(aes(x=mean_p_inf, y=mean_p_vacc, colour=age, shape=race), size=4) +
  geom_linerange(aes(x=mean_p_inf, y=mean_p_vacc, xmin=`2.5%_p_inf`, xmax=`97.5%_p_inf`, colour=age), size=1) +
  geom_linerange(aes(x=mean_p_inf, y=mean_p_vacc, ymin=`2.5%_p_vacc`, ymax=`97.5%_p_vacc`, colour=age), size=1) +
  theme_bw(base_size=14) +
  scale_shape_manual(values=c(15, 16, 17, 21)) +
  xlab("p_inf") + ylab("p_vacc") +
  xlim(0,0.5) + ylim(0,0.5)

## Pull out values to tabulate - updated
summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(
    str_starts(parameter, "p_vacc_given_inf") |
    str_starts(parameter, "p_vacc_given_uninf") |
    str_starts(parameter, "p_vacc") |
    str_starts(parameter, "p_inf") |
    str_starts(parameter, "ratio_")) %>%
  filter(!str_detect(parameter, "p_vacc_prior")) %>%
  mutate(demog_group = as.numeric(str_sub(parameter, -2, -2))) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  as.data.frame() -> to_write_1

#writexl::write_xlsx(to_write_1, path="../two_assay_params_1_v1.xlsx")
writexl::write_xlsx(to_write_1, path="../two_assay_params_1_v2.xlsx")

## Pull out values to tabulate
summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "ratio_")) %>%
  separate(col=parameter, into=c("tmp1","tmp2","tmp3","tmp4","tmp5","demog_group")) %>%
  select(-tmp1, -tmp2, -tmp3, -tmp4, -tmp5) %>%
  mutate(parameter="ratio_p_vacc_p_inf") %>%
  mutate(demog_group = as.numeric(demog_group)) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  as.data.frame()

## -----------------------------------------------------------

##################
## Look at RRRs ##
##################

## Load in model fit
#load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v1.RData")
load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v2.RData")

## Get labels for plotting
dat_for_Stan %>%
  group_by(age_race) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  mutate(label_for_plot = paste0(age_race, " (n=", n, ")")) %>%
  as.data.frame() -> dat_label_age_race

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
  as.data.frame() -> to_write_2

#writexl::write_xlsx(to_write_2, path="../two_assay_params_2_v1.xlsx")
writexl::write_xlsx(to_write_2, path="../two_assay_params_2_v2.xlsx")

## Get labels for plotting
dat_for_Stan %>%
  group_by(age_race) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  mutate(label_for_plot = paste0(age_race, " (n=", n, ")")) %>%
  as.data.frame() -> dat_label_age_race

mcmc_intervals(fit_stan, pars=vars(`p_vacc[1]`:`p_vacc[8]`, `p_inf[1]`:`p_inf[8]`)) +
  ggtitle("Probability of vaccination vs. infection (model with separate p_vacc by infection status & priors)") +
  xlim(0,1)

## RRR
mcmc_intervals(fit_stan, pars=vars("RRR_vs_white[1]","RRR_vs_white[2]","RRR_vs_white[3]","RRR_vs_white[4]"), prob_outer=0.95) +
  scale_x_log10(breaks=scales::breaks_log(6)) +
  scale_y_discrete(labels=dat_label_age_race$age_race[1:4]) +
  # scale_y_discrete(labels=dat_label_age_race$label_for_plot[1:4]) +
  ggtitle("Relative risk ratio of vacc vs. inf\n(ref group: 18-64 & White)") +
  theme_bw(base_size=16) +
  geom_vline(xintercept=1, colour="red") -> p_group_RRR_age1

mcmc_intervals(fit_stan, pars=vars("RRR_vs_white[5]","RRR_vs_white[6]","RRR_vs_white[7]","RRR_vs_white[8]"), prob_outer=0.95) +
  scale_x_log10(breaks=scales::breaks_log(6)) +
  scale_y_discrete(labels=dat_label_age_race$age_race[5:8]) +
  # scale_y_discrete(labels=dat_label_age_race$label_for_plot[5:8]) +
  ggtitle("Relative risk ratio of vacc vs. inf\n(ref group: 65+ & White)") +
  theme_bw(base_size=16) +
  geom_vline(xintercept=1, colour="red") -> p_group_RRR_age2

p_group_RRR_age1 + p_group_RRR_age2

## Save labels
save(dat_label_age_race, file="~/../Box/PHI_scaleit/SCALEIT_2.0/labels.RData")

## -----------------------------------------------------------

## Look at the p_inf_vs_white & p_vacc_vs_white

## Load in model fit
load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v3.RData")

## Get labels for plotting
dat_for_Stan %>%
  group_by(age_race) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  mutate(label_for_plot = paste0(age_race, " (n=", n, ")")) %>%
  as.data.frame() -> dat_label_age_race

## Pull out values to tabulate
summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "p_inf_vs_white")) %>%
  separate(col=parameter, into=c("tmp1","tmp2","tmp3","tmp4","demog_group")) %>%
  select(-tmp1, -tmp2, -tmp3, -tmp4) %>%
  mutate(parameter="p_inf_vs_white") %>%
  mutate(demog_group = as.numeric(demog_group)) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  as.data.frame() -> to_write_3a

summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "p_vacc_vs_white")) %>%
  separate(col=parameter, into=c("tmp1","tmp2","tmp3","tmp4","demog_group")) %>%
  select(-tmp1, -tmp2, -tmp3, -tmp4) %>%
  mutate(parameter="p_vacc_vs_white") %>%
  mutate(demog_group = as.numeric(demog_group)) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  as.data.frame() -> to_write_3b

bind_rows(to_write_3a, to_write_3b) %>%
  select(race, age, parameter, mean, `50%`, `2.5%`, `97.5%`) -> to_write_3

writexl::write_xlsx(to_write_3, path="../two_assay_params_p_inf_p_vacc_relative_to_white.xlsx")

## -----------------------------------------------------------
## -----------------------------------------------------------

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

##################################################
## Table: collate bivariate parameter estimates ##
##################################################

load("~/../Box/PHI_scaleit/SCALEIT_2.0/labels.RData")

## Load in model fit
#load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v1.RData")
load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_v2.RData")

## Pull out values to tabulate - updated
summary(fit_stan)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(
    str_starts(parameter, "p_vacc_given_inf") |
      str_starts(parameter, "p_vacc_given_uninf") |
      str_starts(parameter, "p_vacc") |
      str_starts(parameter, "p_inf") |
      str_starts(parameter, "ratio_")) %>%
  filter(!str_detect(parameter, "p_vacc_prior")) %>%
  mutate(demog_group = as.numeric(str_sub(parameter, -2, -2))) %>%
  mutate(parameter = str_sub(parameter, 1, -4)) %>%
  left_join(dat_label_age_race %>% separate(col=age_race, into=c("age","race"), sep="&") %>% mutate(demog_group=1:8), by="demog_group") %>%
  as.data.frame() -> to_write_1

## Load in model fit
#load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v1.RData")
load("~/../Box/PHI_scaleit/SCALEIT_2.0/pathway_to_immunity_with_RRR_vs_white_v2.RData")

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
  as.data.frame() -> to_write_2

bind_rows(to_write_1, to_write_2) %>%
  select(age, race, parameter, `50%`, `2.5%`, `97.5%`) %>%
  mutate(race = ifelse(race==" Asian", "Asian", race)) %>%
  mutate(race = ifelse(race==" White", "White", race)) %>%
  mutate(race = ifelse(race==" Black or African American", "Black", race)) %>%
  mutate(race = ifelse(race==" Hispanic or Latino", "Latinx", race)) %>%
  rename(median = `50%`) %>%
  rename(lb_2.5 = `2.5%`) %>%
  rename(ub_97.5 = `97.5%`) %>%
  rename(Age = age) %>%
  rename(`Race/Ethnicity` = race) %>%
  as.data.frame() -> to_write_all_bivariate_params

writexl::write_xlsx(to_write_all_bivariate_params, path="../two_assay_params_final.xlsx")
