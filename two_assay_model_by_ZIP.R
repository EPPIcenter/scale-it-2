library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(simstudy)
library(readr)
library(lubridate)
library(tidyverse)
library(patchwork)

source("~/../Box/PHI_scaleit/SCALEIT_2.0/scripts/Stan_functions_v1.R")

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
  ## Remove NAs
  filter(!is.na(seropos_Vitros) & !is.na(seropos_Roche)) %>%
  group_by(age_groups) %>%
  tally()


dat_SCALEIT %>%
  group_by(DateCollected) %>%
  tally()

## Separate the raw univariate data by zcta
dat_SCALEIT %>%
  ## Remove NAs
  filter(!is.na(seropos_Vitros) & !is.na(seropos_Roche)) %>%
  filter(!is.na(zcta)) %>%
  group_by(zcta, seropos_Vitros, seropos_Roche) %>%
  tally() %>%
  ungroup() %>%
  ## Complete missing combos for Vitros & Roche by zcta
  complete(seropos_Vitros, seropos_Roche, nesting(zcta), fill=list(n=0)) %>%
  arrange(zcta, desc(seropos_Vitros), desc(seropos_Roche)) %>%
  ## Remove Roche positive, Vitros negative results
  filter((seropos_Vitros==0 & seropos_Roche==1)==FALSE) %>%
  select(zcta, seropos_Vitros, seropos_Roche, n) %>%
  as.data.frame() -> dat_for_Stan

dat_for_Stan %>%
  group_by(zcta) %>%
  summarise(n=sum(n)) %>%
  arrange(n)

## -----------------------------------------------------------

########################################
## Fit the model - without prior data ##
########################################

fit_stan_1pvacc <- stan(
  model_code=model_stan_Ngroup_1pvacc,
  data=list(
    n_groups=length(unique(dat_for_Stan$zcta)), ## Which one?
    A_plus_C=dat_for_Stan %>% filter(seropos_Vitros==1 & seropos_Roche==1) %>% select(n) %>% pull(),
    B=dat_for_Stan %>% filter(seropos_Vitros==1 & seropos_Roche==0) %>% select(n) %>% pull(),
    D=dat_for_Stan %>% filter(seropos_Vitros==0 & seropos_Roche==0) %>% select(n) %>% pull()))

#save(dat_for_Stan, fit_stan_1pvacc, file="~/../Box/PHI_scaleit/SCALEIT_2.0/bivariate_Stan_results_by_ZIP.RData")

## Load in Stan output
load("~/../Box/PHI_scaleit/SCALEIT_2.0/bivariate_Stan_results_by_ZIP.RData")

## Look at 2D
summary(fit_stan_1pvacc)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "p_inf") | str_starts(parameter, "p_vacc\\[")) %>%
  separate(col=parameter, into=c("tmp","parameter","zcta")) %>%
  select(-tmp) %>%
  mutate(parameter=paste0("p_", parameter)) %>%
  mutate(zcta = rep(unique(dat_for_Stan$zcta), times=2)) %>%
  as.data.frame() %>%
  pivot_wider(names_from=parameter, values_from=mean:Rhat) %>%
  mutate(n_samples = dat_for_Stan %>% group_by(zcta) %>% summarise(n=sum(n)) %>% select(n) %>% pull()) %>%
  filter(n_samples >= 10) %>%
  ggplot() +
  geom_linerange(aes(x=mean_p_inf, y=mean_p_vacc, xmin=`2.5%_p_inf`, xmax=`97.5%_p_inf`, colour=factor(zcta)), size=0.05) +
  geom_linerange(aes(x=mean_p_inf, y=mean_p_vacc, ymin=`2.5%_p_vacc`, ymax=`97.5%_p_vacc`, colour=factor(zcta)), size=0.05) +
  geom_abline(linetype="dotted") +
  geom_point(aes(x=mean_p_inf, y=mean_p_vacc, colour=factor(zcta), size=n_samples)) +
  #ggrepel::geom_text_repel(aes(x=mean_p_inf, y=mean_p_vacc, label=zcta), min.segment.length=0, direction="both", force=2) +
  theme_bw(base_size=14) +
  xlab("Probability of infection") + ylab("Probability of vaccination")

#######################
## Put them on a map ##
#######################

library(ggplot2)
library(ggmap)
library(ggpubr)
library(maps)
library(mapproj)
library(mapdata)
library(sf)
library(tidyverse)
library(geojsonio)
library(broom)
library(viridis)

spdf <- geojson_read("~/../Box/PHI_scaleit/SCALEIT_2.0/Bay Area ZIP Codes.geojson",  what = "sp")
spdf<-spdf[substr(spdf$zip, 1,3)  %in% c("941") , ]
spdf<-spdf[spdf$zip != "94128" , ]


spdf_fortified <- tidy(spdf, region = "zip")
spdf_fortified$id<-as.numeric(spdf_fortified$id)

## Clean modelled output
summary(fit_stan_1pvacc)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  filter(str_starts(parameter, "p_inf") | str_starts(parameter, "p_vacc\\[")) %>%
  separate(col=parameter, into=c("tmp","parameter","zcta")) %>%
  select(-tmp) %>%
  mutate(parameter=paste0("p_", parameter)) %>%
  mutate(zcta = rep(unique(dat_for_Stan$zcta), times=2)) %>%
  as.data.frame() %>%
  pivot_wider(names_from=parameter, values_from=mean:Rhat) %>%
  mutate(n_samples = dat_for_Stan %>% group_by(zcta) %>% summarise(n=sum(n)) %>% select(n) %>% pull()) %>%
  mutate(`50%_p_inf` = ifelse(n_samples < 10, NA, `50%_p_inf`)) %>%
  mutate(`50%_p_vacc` = ifelse(n_samples < 10, NA, `50%_p_vacc`)) %>%
  as.data.frame() -> zcta_modelled

spdf_fortified = spdf_fortified %>%
  left_join(. ,  as.data.frame(zcta_modelled), by=c("id"="zcta"))

ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill = `50%_p_inf`, x = long, y = lat, group = group)) +
  theme_void() +
  coord_map()+
  #theme(legend.position = "none")+
  labs_pubr()+
  rremove("axis.text")+
  rremove("xlab")+
  rremove("ylab")+
  scale_fill_distiller(name="value", palette ="RdBu")+
  ggtitle("a. Probability of prior infection") -> plot_inf

ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill = `50%_p_vacc`, x = long, y = lat, group = group)) +
  theme_void() +
  coord_map()+
  #theme(legend.position = "none")+
  labs_pubr()+
  rremove("axis.text")+
  rremove("xlab")+
  rremove("ylab")+
  scale_fill_distiller(name="value", palette = "BrBG")+
  ggtitle("b. Probability of vaccination") -> plot_vacc


ggarrange(plot_inf, plot_vacc,
          #labels = c("a", "b"),
          ncol = 2, nrow = 1, align="hv")

ggsave("~/../Box/PHI_scaleit/SCALEIT_2.0/modeled_map.tiff", units="in", width=9, height=4, dpi=300, compression = 'lzw')
ggsave("~/../Box/PHI_scaleit/SCALEIT_2.0/modeled_map.png", units="in", width=9, height=4, dpi=300)


#writexl::write_xlsx(as.data.frame(zcta_modelled) %>% filter(n_samples >= 10), path="~/../Box/PHI_scaleit/SCALEIT_2.0/bivariate_ZIP_estimates.xlsx")
